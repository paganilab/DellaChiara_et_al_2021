#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                  RNA-SEQ ORGANOIDS   
========================================================================================
 @#### Authors
 Raoul Bonnal <raoul.bonnal@ifom.eu>
 Michaela Fakiola <michaela.fakiola@ifom.eu>
 Federica Gervasoni <federica.gervasoni@ifom.eu>
----------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------
Pipeline overview:
 - STEP 1 - FastQC: quality control of sequencing data
 - STEP2 - Removal of adapters and trimming of reads using bbduk
 - STEP 3 - FastQC: Quality control assessment of fastqc data after trimming
 - STEP4 - Alignment of reads using STAR and indexing using Sambamba
 - STEP5 - Generate normalised coverage track file using bamCoverage.
 - STEP6 - Quantification of reads
----------------------------------------------------------------------------------------
*/

/* Example of command line call

nextflow run main.nf --csv sampleteset.csv -w /bulk_rnaseq/work -resume

*/

HTseqGenes = file(params.genome.HTseqGenes)
adapters = params.adapters 
intronBed = file(params.genomesdir + '/' + params.intron_bed.local.path + '/' + params.intron_bed.local.name)
gtf = file(params.genomesdir + '/' + params.gtf.local.path + '/' + params.gtf.local.name)
rRNA = file(params.genomesdir + '/' + params.rRNA.local.path + '/' + params.rRNA.local.name)
mtRNA = file(params.genomesdir + '/' + params.mtRNA.local.path + '/' + params.mtRNA.local.name)
index_star=file(params.genomesdir + '/' + params.star_index)

readsToAggregate = Channel.create()

readsToProcess = Channel.create()

Channel
    .from(file(params.csv))
    .splitCsv(header: true)
    .map { sample -> 
    [sample.origin, sample.sample, file(sample.r1), file(sample.r2)] }
    .groupTuple(by:[0,1])
    .choice(readsToAggregate, readsToProcess) { sample ->
    sample[2].size() > 1 ? 0 : 1 }

process aggregateReads {
    cpus 2
    memory '4GB'
    input:
        set val(origin), val(sample), file(r1), file(r2) from readsToAggregate

    output:
        set val(origin), val(sample), file("${sample}_R1.fastq.gz"), file("${sample}_R2.fastq.gz") into readsAggregated

    script:
    """
    zcat ${r1} | gzip  > ${sample}_R1.fastq.gz
    zcat ${r2} | gzip  > ${sample}_R2.fastq.gz
    """
}

readsToProcess
    .mix(readsAggregated)
    .into { raw_reads_fastqc; raw_reads_trimming }

process fastqc {

/*
 * STEP 1 - FastQC: quality control of sequencing data
 */
    tag "fastqc-${sample}"
    
    cpus 3
    memory '4GB'
    publishDir "${params.outdir}/${sample}/fastqc", mode: 'copy'

    input:
	set val(origin), val(sample), file(r1), file(r2)  from raw_reads_fastqc

    output:
	file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    """
}

process trimming { 
/*
 * STEP2 - Removal of adapters and trimming of reads using bbduk
 */
    tag "Trimming-${sample}"
    cpus 8
    publishDir "${params.outdir}/${sample}/trimming", mode: 'copy'

    input:
	set val(origin), val(sample), file(r1), file(r2) from raw_reads_trimming

    output:
	set val(sample), file("${sample}_R1_trim.fastq.gz"), file("${sample}_R2_trim.fastq.gz") into reads_trimmed, reads_trimmed_qc

    script:
    """
    bbduk.sh in1=${r1} in2=${r2} out1=tmp1.fastq.gz out2=tmp2.fastq.gz ref=${adapters} ktrim=r k=23 mink=11 hdist=1 t=${task.cpus} tpe tbo qin=33
    bbduk.sh in1=tmp1.fastq.gz in2=tmp2.fastq.gz out1=${sample}_R1_trim.fastq.gz out2=${sample}_R2_trim.fastq.gz qtrim=rl trimq=20 minlen=50 t=${task.cpus} qin=33
    """
}

process fastqcPostTrimming {
/*
 * STEP 3 - FastQC: Quality control assessment of fastqc data after trimming
 */
    tag "fastqcTrim-${sample}"
    cpus 3
    memory '4GB'
    publishDir "${params.outdir}/${sample}/fastqcTrimming", mode: 'copy'

    input:
	set val(sample), file(r1), file(r2) from reads_trimmed_qc

    output:
	file '*_fastqc.{zip,html}' into fastqc_trimmed_results

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    """
}

process starMapping {
/*
 * STEP4 - Alignment of reads using STAR and indexing using Sambamba
 */
  tag "StarMapping-${sample}"
  cpus 20
  memory '40G'

  publishDir "${params.outdir}/${sample}/mapped", mode: 'copy'

  input:
    file(index_star)
    set val(sample), file(r1), file(r2) from reads_trimmed
  output:
 	set val(sample), file("${sample}_Aligned.sortedByCoord.out.bam"), file("${sample}_Aligned.sortedByCoord.out.bam.bai") into bams_to_quantification, bams_to_wig, bams_to_mtrna, bams_to_rrna
    set val(sample), file("${sample}_Aligned.sortedByCoord.out.bam"), file("${sample}_Aligned.sortedByCoord.out.bam.bai"), file("${sample}_Aligned.sortedByCoord.out.flagstats") into bams_to_retained_introns    
    set val(sample), file("${sample}_Log.final.out") into mappedLogFinal
    set val(sample), file("${sample}_Log.out") into mappedLog
    set val(sample), file("${sample}_Log.progress.out") into mappedProgress
    set val(sample), file("${sample}_ReadsPerGene.out.tab") into  mappedReadsPerGene
    set val(sample), file("${sample}_SJ.out.tab") into mappedSJ
    set val(sample), file("${sample}_Unmapped.*") into mappedUnmapped
    set val(sample), file("${sample}_Aligned.sortedByCoord.out.*stats") into mappedStats
    
    script:
    """
STAR --genomeDir ${index_star} \
     --runThreadN ${task.cpus} \
     --readFilesIn ${r1} ${r2} \
     --readFilesCommand zcat \
     --genomeLoad LoadAndRemove \
     --outFileNamePrefix ${sample}_ \
     --outReadsUnmapped Fastx \
     --outSAMtype BAM SortedByCoordinate \
     --limitBAMsortRAM  ${params.samtools.memory} \
     --alignIntronMax 1000000 \
     --quantMode GeneCounts \
     --outFilterMismatchNmax 9 \
     --outFilterMultimapNmax 20 \
     --alignSJoverhangMin 8 \
     --alignSJDBoverhangMin 1 \
     --alignMatesGapMax 1000000

   sambamba index ${sample}_Aligned.sortedByCoord.out.bam -t ${task.cpus}
   samtools flagstat ${sample}_Aligned.sortedByCoord.out.bam > ${sample}_Aligned.sortedByCoord.out.flagstats
   samtools idxstats ${sample}_Aligned.sortedByCoord.out.bam > ${sample}_Aligned.sortedByCoord.out.idxstats
   """
}

process bigWig{
/*
 * STEP5 - Generate normalised coverage track file using bamCoverage. 
 * Multimapped reads are filtered so that the tracks more closely match the 
 * read quantification data. Separate tracks for forward and reverse transcripts
 * are generated. Two approaches were followed: using samtools to pre-filter or
 * directly in bamCoverage. Since the output files from the two approaches are 
 * very similar, the samtools approach is masked. 
 */
  tag "BigWig-${sample}"
  cpus 8
  memory '40G'

  publishDir "${params.outdir}/${sample}/bigwig", mode: 'copy'

  input:
    set val(sample), file(bam), file(bai) from bams_to_wig
  output:
    set val(sample), file("${sample}_forward.bw"), file("${sample}_reverse.bw") into bigwigs

  script:
  """
  bamCoverage -b ${bam} -o ${sample}_forward.bw --filterRNAstrand forward --normalizeTo1x 3049315783 --minMappingQuality 10 -p ${task.cpus}
  bamCoverage -b ${bam} -o ${sample}_reverse.bw --filterRNAstrand reverse --normalizeTo1x 3049315783 --minMappingQuality 10 -p ${task.cpus}
  """
}

process quantification {

/*
 * STEP6 - Quantification of reads
 */    

    cpus 8
    publishDir "${params.outdir}/${sample}/quantification", mode: 'copy'
    input:
	file(gtf)
	set val(sample), file("data.bam"), file("data.bam.bai") from bams_to_quantification
	
    output:
	set val(sample), file("${sample}_counts.txt") into quantifications, mappedReadsPerGeneToAggegate 
        set val(sample), file("${sample}_counts.txt.summary") into quantification_stats

    script:
    """
    featureCounts -p -B -C -t exon -g gene_id -a ${gtf} -o ${sample}_counts.txt -s 2 -T ${task.cpus} data.bam 
    """
}

HTseqCounts = mappedReadsPerGeneToAggegate.map{ sample, count -> count}

process aggregateHTseqCounts {
    cpus 2
    memory '2G'
    publishDir "${params.outdir}/", mode: 'copy'
    input:
    file(counts) from HTseqCounts.collect()
    file(HTseqGenes)

   output:
    file('HTseqCounts.tsv')

    script:
    """
    echo "ENSG" > HTseqCounts.tsv
    cat ${HTseqGenes} >> HTseqCounts.tsv
    for f in *_counts.txt
    do
      sample=\$(basename \${f} _counts.txt)
      echo \${sample} > ftmp
      tail -n +3 \${f} |  awk '{print \$7}' >> ftmp
      paste HTseqCounts.tsv ftmp > HTseqCounts.tmp
      mv HTseqCounts.tmp HTseqCounts.tsv 
    done
    rm ftmp
    """
}

process count_mtrna {
/*
 * Count mitochondrial rna
 */
    tag "countMtRNA-${sample}"
    cpus 8
    publishDir "${params.outdir}/${sample}/quantification", mode: 'copy'

    input:
	file(mtRNA)
        set val(sample), file("data.bam"), file("data.bam.bai") from bams_to_mtrna

    output:
	set val(sample), file("${sample}_mtRNAsplit*") into mtrnas

    script:
    """
    split_bam.py -i data.bam -r ${mtRNA} -o ${sample}_mtRNAsplit
    """
}

process count_rrna {
/*
 * count ribosomal rna
 */
    tag "countRRNA-${sample}"
    publishDir "${params.outdir}/${sample}/quantification", mode: 'copy'
    cpus 8

    input:
	set val(sample), file("data.bam"), file("data.bam.bai") from bams_to_rrna
    file(rRNA)

    output:
	set val(sample), file("${sample}_rRNAsplit*") into rrnas	
    
    script:
    """
    split_bam.py -i data.bam -r ${rRNA} -o ${sample}_rRNAsplit
    """
}
  
process  count_retained_introns {
/*
 *  count retained intron events
 */
    
    tag "countRetainedIntrons-${sample}"
    cpus 8
    publishDir "${params.outdir}/${sample}/quantification", mode: 'copy'
    input:
	file(intronBed)
    set val(sample), file("data.bam"), file("data.bam.bai"), file("data.flagstats") from bams_to_retained_introns

    output:
	file("${sample}/*") into countsRetainedIntrons

    script:
    """
    mappedreads=\$(sed -n '5,5p' data.flagstats | cut -f1 -d ' ')
    iread.py data.bam ${intronBed} -t \${mappedreads} -o ${sample}
    """
}
