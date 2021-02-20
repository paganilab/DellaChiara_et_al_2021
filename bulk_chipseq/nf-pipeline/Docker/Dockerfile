FROM pipelines/base:0.4

LABEL org.bioinformatics.authors="Federica Gervasoni, Raoul Jean Pierre Bonnal"
LABEL maintainer="federica.gervasoni@ifom.eu, raoul.bonnal@ifom.eu"
LABEL description="Docker image containing all requirements for ChIP-seq pipeline"
LABEL version="v0.1.5.1b0.4"

RUN apt-get update && apt-get install -y liblapack3 libboost1.62-dev gawk

RUN conda install -y \
    bedtools=2.27.1 \
    bowtie=1.2.2 \
    bowtie2=2.3.4.1 \
    bioconductor-biomart=2.34.2 \
    bioconductor-bsgenome=1.46.0 \
    bioconductor-bsgenome.hsapiens.ucsc.hg19=1.4.0 \
    bioconductor-bsgenome.mmusculus.ucsc.mm10=1.4.0 \
    bioconductor-chippeakanno=3.12.0 \
    bioconductor-ensembldb=2.2.2 \
    bioconductor-genomicranges=1.30.3 \
    bioconductor-genomicfeatures=1.30.3 \
    bioconductor-org.mm.eg.db=3.5.0 \
    bioconductor-rsamtools=1.30.0 \
    bioconductor-rtracklayer=1.38.3 \
    bioconductor-shortread=1.36.0 \
    fastqc=0.11.7 \
    multiqc=1.5 \
    samtools=1.8 \
    picard=2.18.6 \
    r-spp=1.15.2 \
    r-snowfall=1.84-6.1 \
    r-catools=1.17.1 \
    trimmomatic=0.36 \
    ucsc-bedgraphtobigwig=366 \
    sra-tools=2.9.1 &&\
    conda clean --all

RUN conda install -y -c conda-forge openjdk=8.0.152 &&\
    conda install -y -c anaconda fontconfig=2.13.1 \
                                 font-ttf-dejavu-sans-mono=2.37 &&\
    conda clean --all

WORKDIR /opt
RUN git clone https://github.com/kundajelab/phantompeakqualtools &&\
    cd phantompeakqualtools &&\
    git checkout 916b09f1b7edba4d47d61ce4930c78fd47e807d9 &&\
    sed -i '1s/^/#\!\/usr\/bin\/env Rscript\nlibrary(caTools)\n/' run_spp.R &&\
    chmod 755 run_spp.R &&\
    Rscript --slave --no-save --no-restore-history -e 'install.packages("spp_1.14.tar.gz",repos="https://cloud.r-project.org/")'

RUN conda create -n macs2 python=2.7 &&\
    /bin/bash -c "conda activate macs2; conda install macs2=2.1.0;conda clean --all; conda deactivate macs2"

RUN conda create -n deeptools &&\
    /bin/bash -c "conda activate deeptools; conda install deeptools=3.0.2;conda clean --all; conda deactivate deeptools"

WORKDIR /opt

RUN curl -L -o BBMap_38.16.tar.gz "https://sourceforge.net/projects/bbmap/files/BBMap_38.16.tar.gz/download" &&\
    tar xzvf BBMap_38.16.tar.gz &&\
    rm BBMap_38.16.tar.gz

ENV PATH=/opt/bbmap:/opt/phantompeakqualtools:${PATH}
