//pipeline: RNA-seq
/*

Running the PRE script for setting up all files required by the regular RNA-seq pipeline.


Example command line:
nextflow run pre-main.nf -c nextflow.config -w work

*/

Channel.from(
    ['',params.rRNA],
    ['',params.mtRNA],
    ['',params.intron_bed],
    ["genome:HG38v25", params.fasta],
    ["genome:HG38v25",params.gtf])
    .set { files }

process getFile {
    tag "${name}"
    publishDir "${params.genomesdir}/${path}", overwrite: true, mode: 'copy'
    input:
	set val(tag), val(data) from files
    output:
        set val(tag), file("*") into ofiles
    script:
        remote = data.remote.url
        path = data.local.path
        name = data.local.name
    format = data.remote.format
    if (format=='text') {
        """
wget -O ${name} "${remote}"
"""
    }
    else if (format=='gzip') {
	    """
wget -O tmp "${remote}"
zcat tmp > ${name}
rm tmp
"""
    }
    else {
	    """
wget -O ${name} "${remote}"
"""
    }
}

ofiles.filter({ it[0]=~/genome:.*/ })
    .groupTuple()
    .map {
	  isFirstFasta = it[1][0].getExtension()=='fasta'
	  [it[0],
	   isFirstFasta ? it[1][0] : it[1][1],
	   isFirstFasta ? it[1][1] : it[1][0]]}
    .set {genome}
  
process index {
    tag "Create_Index_${reference}"

    cpus 20
    memory '40G'

    publishDir "${params.genomesdir}/${params.index_star}", mode: 'move'

    input:
	set val(reference), file(fasta), file(gtf) from genome

    output:
	set val(reference), file('*') into genomeIndexes

    script:
            """
STAR --runThreadN ${task.cpus} \
     --runMode genomeGenerate \
     --genomeDir . \
     --sjdbGTFfile ${gtf} \
     --genomeFastaFiles ${fasta} \
     --sjdbOverhang 74
"""
}
