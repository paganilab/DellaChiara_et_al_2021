#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                  CHIP-SEQ ORGANOIDS   
========================================================================================
 [organoids]ChIP-seq Pipeline. Started Jan 2018.
 @#### Authors
 Federica Gervasoni <federica.gervasoni@ifom.eu>
 Raoul JP Bonnal <raoul.bonnal@ifom.eu>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - STEP 1: FastQC: quality control of sequencing data
 - STEP 2: align with bowtie, remove multimapping and states
 - STEP 3: Remove duplicates and stats
 - STEP 4: Stats_remove_duplicates
 - STEP 5: Phantompeakqualtools
 - STEP 6: MACS2 peak calling
 - STEP 7: Remove blacklisted regions and scaffolds
 - STEP 8: bw creation
 - STEP 9: deepTools
 - STEP 10: MultiQC 
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    =========================================
     ChIPseq : ChIP-Seq Best Practice v${version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --genome hg38 --nomodel --reads fastqfiles_test.csv --macsconfig configfile_test.csv -w nfwd -resume

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes).
      --genome                      Name of iGenomes reference
      --macsconfig                  Configuration file for peaking calling using MACS. Format: ChIPSampleID,CtrlSampleID,AnalysisID
      --profile                     Hardware config to use: Docker

    Options:
      --broad                       Run MACS with the --broad flag
      --nomodel                     Run MACS with the --nomodel --extsize 200
      --sra                         If reference data must be downloaded from SRA/NCBI

    Trimming options
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.

    Other options:
      --outdir                      The output directory where the results will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Pipeline version
version = '1'

params.genomes['hg38'].bowtie_index = "${params.genomesbase}/Index"
params.genomes['hg38'].blacklist = file("${params.genomesbase}/blacklists/hg38.blacklist.bed")
params.genomes['hg38'].adapters= file("${params.genomesbase}/adapters.fa") 
params.genomes['hg38'].gtf = file("${params.genomesbase}/gencode.v25.basic.annotation.gtf")
params.genomes['hg38'].fasta = file("${params.genomesbase}/fasta/GRCh38.primary_assembly.genome.fa")
params.genomes['hg38'].chrome_size = "${params.genomesbase}/hg38.chrom.sizes_macs2"

println(params.genomes['hg38'])

// Configurable variables
params.name = false
params.project = false
params.genome = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bowtie_index = params.genome ? params.genomes[ params.genome ].bowtie_index ?: false : false
params.bowtie_name = params.genome ? params.genomes[ params.genome ].bowtie_name ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.chrome_size = params.genome ? params.genomes[ params.genome ].chrome_size ?: false : false
params.genome_size = params.genomes[ params.genome ].size
params.adapters = params.genome ? params.genomes[ params.genome ].adapters ?: false : false
params.reads ="fastqfiles_organoids_bulk_chipseq_test.txt"
params.notrim = false
params.broad = false
params.nomodel = params.nomodel ? params.nomodel : false  
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
params.outdir = params.outdir ? params.outdir : './results'
params.singleEnd="singleEnd"

blacklist = file(params.blacklist)
    if( !blacklist.exists() ) exit 1, "Blacklist file not found: ${params.blacklist}"

gtf = file(params.gtf)

bowtie_index = file(params.bowtie_index)

reads_file = file(params.reads)

chrome_size = file(params.chrome_size)

def data = [:]
def design = [:]

Channel.fromPath(params.macsconfig)
    .splitCsv(header: true)
    .subscribe { row ->

    design.put(row.sample,["sample_id": row.sample, "ctrl_id":  row.input, "modification": row.modification, "type": row.chiptype.toLowerCase() ])
}

def get_modification(experimentalDesign, sampleName){    
    experimentalDesign.get(sampleName)['modification']
}

def get_modification_path(experimentalDesign,sampleName){
    
    if (experimentalDesign.containsKey(sampleName)){
    
	get_modification(experimentalDesign,sampleName) + "/" + sampleName
    } else {
	"input/"+sampleName
    }
}

def genome_size = false
if (params.genome_size ){
    genome_size = params.genome_size
} else {
    log.warn "No reference genome size, use the default 3049315783"
    genome_size = 3049315783
}

// Header log info
log.info "========================================="
log.info " ChIP-Seq Best Practice v${version}"
log.info "========================================="
def summary = [:]
summary['Reads']               = params.reads
summary['Data Type']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']              = params.genome
if(params.bowtie_index)  summary['Bowtie Index'] = params.bowtie_index
else if(params.fasta) summary['Fasta Ref'] = params.fasta
if(params.gtf)  summary['GTF File'] = params.gtf
summary['Genome Chrome Size']       = params.chrome_size
summary['MACS Config']         = params.macsconfig
summary['MACS nomodel extsize'] = params.nomodel
summary['Blacklist filtering'] = params.blacklist
if( params.blacklist_filtering ) summary['Blacklist BED'] = params.blacklist
summary['Extend Reads']        = "Default 200 bp"
summary['Current home']        = "$HOME"
summary['Current user']        = "$USER"
summary['Current path']        = "$PWD"
summary['Working dir']         = workflow.workDir
summary['Output dir']          = params.outdir
summary['Script dir']          = workflow.projectDir

if( params.notrim ){
    summary['Trimming Step'] = 'Skipped'
}
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "===================================="

// Check that Nextflow version is up to date enough
nf_required_version = '0.25.0'
try {
    if( ! nextflow.version.matches(">= $nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}

params.filter = ".*"

readsToAggregate = Channel.create()

readsToProcess = Channel.create()

Channel
    .fromPath(params.reads)
	.splitCsv(header: true)
    .map { sample ->
    def name = sample.name
    def fullpath = sample.path
    if ( workflow.profile == 'aws' ) {
	fullpath = sample.awspath
    }
    [name, file(fullpath)] }
    .filter{ (it[0] ==~ /${params.filter}/) || (it[1] ==~ /${params.filter}/) }
    .groupTuple()
    .choice(readsToAggregate, readsToProcess) { sample ->
    sample[1].size() > 1 ? 0 : 1 }

process aggregateReads {
    cpus 2
    memory '4GB'
    input:
	set val(name), file(reads) from readsToAggregate

    output:
	set val(name), file("${name}.fastq.gz") into readsAggregated

    script:
    """
    zcat ${reads} | gzip  > ${name}.fastq.gz
    """	
}

readsToProcess
    .mix(readsAggregated)
    .into { raw_reads_fastqc; raw_reads_bowtie}

/*
 * STEP 1 - FastQC: quality control of sequencing data
 */
process fastqc {
    cpus 3
    memory '4GB'
    publishDir "${params.outdir}/${path}/fastqc", mode: 'copy'

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file '*_fastqc.{zip,html}' into fastqc_results, in_fastqc_results

    script:
    path=get_modification_path(design,name)
    """
    fastqc -t ${task.cpus} $reads 
    """
}

fastqc_results.subscribe { println "Received: " + it }

/*
 * STEP 2 - align with bowtie, remove multimapping and states
 */
process bowtie {
    cpus 8
    memory '10 GB'
    publishDir "${params.outdir}/${path}/bowtie", mode: 'copy', pattern: '*.bowtie.stats'
    input:
	file index from bowtie_index
	set val(name), file(reads) from raw_reads_bowtie
	
    output:
	set val("${name}"),file("${name}.bam"), file("${name}.bowtie.stats") into in_remove_duplicates, in_remove_duplicates_MQ
    script:
	def samtools_cpus = 1
        def bowtie_cpus_keep = 4
	path=get_modification_path(design,name)
    """
    zcat ${reads} | bowtie -S --threads ${task.cpus-bowtie_cpus_keep} -m 1 --best --strata -v 3 ${index.name}/${params.bowtie_name} - 2> ${name}.bowtie.stats | samtools view -bS - | samtools sort -n -@ ${samtools_cpus} -m ${params.samtools.memory} - -T ${name} -O BAM -o ${name}.bam
    """
}

/*
 * STEP 3 - Remove duplicates and stats 
 */
process remove_duplicates {
    cpus 6
    memory '12 GB'
    publishDir path: "${params.outdir}/${path}/remove_duplicates", mode: 'copy'
    input:
	set val(name),file(bam), file(stderr) from in_remove_duplicates

    output:
        set val("${name}"),file("${name}_no_dup.bam"),file("${name}_no_dup.bam.bai") into in_phantompeak, in_deeptools, in_macs_without_duplicates, in_rm_dup_stats

    script:
	path=get_modification_path(design,name)
	
    """
   # starts from reads sorted by name and does: Fill in mate coordinates, ISIZE and mate related flags from a name-sorted alignment.
   # Then we resort by coordinates to proceed with the markdup step, in which the duplicates are marked and removed and a report statistic is created.
   samtools fixmate -r -m ${bam} -| samtools sort -m ${params.samtools.memory} -@ 4 - | samtools markdup -r -s - ${name}_no_dup.bam
   samtools index ${name}_no_dup.bam
    """
}

in_macs_modifications = in_macs_without_duplicates.map{
    data.put(it[0], it)
    it[0]}
    .filter {design.containsKey(it)}

process buildMacsComparison {
    executor 'local'
    memory '2 GB'
    cpus 1
    cache false
    input:
	each sample_id from in_macs_modifications.collect()

    output:
	val(chipCMP) into chipCMPs
 
    exec:
	def sample_info = design.get(sample_id)
	def modification = data.get(sample_id)
        def ctrl = null
        def ctrl_bam=null
        def ctrl_bai=null
        def ctrl_id = sample_info.get("ctrl_id")
        if (ctrl_id!="no"){
           ctrl = data.get(ctrl_id)
           ctrl_bam = ctrl[1]
           ctrl_bai = ctrl[2]
        }
           
    chipCMP =[sample_id, ctrl_id, sample_info.get("type"), modification[1], modification[2], ctrl_bam, ctrl_bai]
    println(chipCMP)
}

/*
 * STEP 4 - Stats_remove_duplicates
 */
process statsRemoveDuplicates {
    cpus 1
    memory '6 GB'
    publishDir "${params.outdir}/${path}/stats_remove_duplicates", mode: 'copy'
    input:
        set val(name),file(bam),file(bai) from in_rm_dup_stats

    output:
       file '*.txt' into rm_dup_stats, rm_dup_stats_MQ

    script:
       path=get_modification_path(design,name)

    """
    samtools flagstat ${bam} > ${name}_flagstat.txt
    """
}


/*
 * STEP 5 - Phantompeakqualtools
 */
process phantompeakqualtools {
    cpus 5
    memory '6 GB'
    publishDir "${params.outdir}/${path}/phantompeakqualtools", mode: 'copy'
    input:
	set val(name),file(bam),file(bai) from in_phantompeak

    output:
       file '*.pdf' into cross_correlation_pdf
       file '*_crossCorrQC.txt' into cross_correlation_qc

    script:
       path=get_modification_path(design,name)
	
    """
    run_spp.R -c=${bam} -s=-200:5:1000 -savp -out=${name}_crossCorrQC.txt
    """
}

/*
 * STEP 6 - MACS2 peak calling
 */
process macs {
    publishDir "${params.outdir}/${path}/macs", mode: 'copy', pattern: '*.{r,xls,stats}'
    memory '8 GB'
    input:
	    set val(sample_name), val(ctrl_name), val(chipType),file(sample_bam),file(sample_bai),file(ctrl_bam),file(ctrl_bai) from chipCMPs 

    output:
       set val(sample_name),file('*.{r,bed,xls,stats}') into macs_results
    set val(sample_name), val(chipType), file("${sample_name}_peaks.${fileext}"),stdout into out_macs_peak
       set val(sample_name),file('*.bdg') into out_macs_bdg
       
    script:
    broad = chipType=="sharp" ? '': "--broad"
    fileext = chipType=="sharp" ? 'narrowPeak' : 'broadPeak'
    nomodel = params.nomodel ? "--nomodel --extsize 200" : ''
    path=get_modification_path(design,sample_name)
    with_input= ctrl_name=="no" ? '' : "-c ${ctrl_bam}"    
    """
    source activate macs2
    macs2 callpeak -t ${sample_bam} ${with_input} ${broad} -f BAM -g ${genome_size} -n ${sample_name} ${nomodel} -B -q 0.01 2> >(tee -a ${sample_name}.macs.stats >&2) 1>${sample_name}.macs.stdout
    wc -l < ${sample_name}_peaks.${fileext}
    """
}

/*
 * STEP 7 - Remove blacklisted regions and scaffolds
 */
process RmBlacklisted {
    publishDir "${params.outdir}/${path}/macs", mode: 'copy'
    memory '2 GB'
    input:
	set val(name), val(chipType), file(peaks), val(nPeaks) from out_macs_peak
    file(bl) from blacklist
  
    output:
	file "${name}_noBL_noUN.${fileext}" into in_bw_creation

    when:
	(nPeaks as Integer) > 0

    script:
      fileext = chipType=='broad' ? "broadPeak" : "narrowPeak"
      def blacklisted = params.blacklist ? "-b ${bl}" : ''
      path=get_modification_path(design,name)
    """
    bedtools intersect -a $peaks $blacklisted -v | grep chr > ${name}_noBL_noUN.${fileext}
    """
}

/*
 * STEP 8 - bw creation
 */
process bwCreation {
    publishDir "${params.outdir}/${path}/macs", mode: 'copy'
    memory '48 GB'
    input:
	set val(name), file(peaks_bdg) from out_macs_bdg
        file chrome_size

    output:
      file "${name}_FE.bw" into out_bw

    script:

    def treat_pileup = "${name}_treat_pileup.bdg"
    def control_lambda = "${name}_control_lambda.bdg"
	path=get_modification_path(design,name)

    """
    source activate macs2
    macs2 bdgcmp -t $treat_pileup -c $control_lambda -o ${name}_FE.bdg -m FE
    source deactivate
    LC_COLLATE=C sort -k1,1 -k2,2n ${name}_FE.bdg > ${name}_FE_sort.bdg
    bedGraphToBigWig ${name}_FE_sort.bdg $chrome_size ${name}_FE.bw 
    """
}

/*
 * STEP 9 - deepTools
 */
process deepTools {
    cpus 4
    memory "24 GB"
    publishDir "${params.outdir}/${path}/deepTools", mode: 'copy'
    input:
        set val(name),file(bam),file(bai) from in_deeptools
        file annotation from gtf
        val ref_deepTools from genome_size
    output:
      file '*.{pdf,gz,bw}' into deepTools

    script:
    path=get_modification_path(design,name)

    """
    source activate deeptools
    export TMPDIR=\$PWD
    bamCoverage --numberOfProcessors ${task.cpus} -v -b $bam -o ${name}_norm.coverage.bw --normalizeUsing RPGC --effectiveGenomeSize $ref_deepTools --extendReads 200 --binSize 1
    computeMatrix reference-point --numberOfProcessors ${task.cpus} --regionsFileName $annotation --scoreFileName ${name}_norm.coverage.bw --outFileName ${name}.TSS.gz --upstream 3000 --downstream 3000
    computeMatrix scale-regions --numberOfProcessors ${task.cpus} --regionsFileName $annotation --scoreFileName ${name}_norm.coverage.bw --outFileName ${name}.genebody.gz --regionBodyLength 6000 --upstream 3000 --downstream 3000
    plotProfile --matrixFile ${name}.TSS.gz --outFileName ${name}.TSS.profile.pdf
    plotProfile --matrixFile ${name}.genebody.gz --outFileName ${name}.genebody.profile.pdf
    plotHeatmap --matrixFile ${name}.TSS.gz --outFileName ${name}.TSS.heatmap.pdf --colorMap seismic --zMin 0 --zMax 10 --heatmapHeight 20
    plotHeatmap --matrixFile ${name}.genebody.gz --outFileName ${name}.genebody.heatmap.pdf --colorMap seismic --zMin 0 --zMax 10 --heatmapHeight 20
    """
}

/*
 * STEP 10 - MultiQC
 */
process multiqc {
    cpus 1
    memory "4 GB"
    publishDir "${params.outdir}/multiQC", mode: 'copy'

    input:
    file ('fastqc/*') from in_fastqc_results.collect()
    file ('bowtie/*') from in_remove_duplicates_MQ.collect()
    file ('stats_remove_duplicates/*') from rm_dup_stats_MQ.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_data

    script:
    """
    multiqc .  
    """
}
