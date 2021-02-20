FROM  bionic-base:0.1

#
# How to build
# docker build -t organoids-rna-seq:0.1 .
# docker tag organoids-rna-seq:0.1 localhost:5000/organoids-rna-seq:0.1
# docker push localhost:5000/organoids-rna-seq:0.1
#

LABEL org.bioinformatics.authors="Michaela Fakiola, Federica Gervasoni, Raoul Jean Pierre Bonnal"
LABEL org.bioinformatics.maintainer="raoul.bonnal@ifom.eu"
LABEL org.bioinformatics.description="Docker image containing all requirements for RNA Sequencing pipeline" 
LABEL software.iread="https://www.biorxiv.org/content/early/2017/10/04/135624www.libpls.net/iread"

ENV DEBIAN_FRONTEND=noninteractive
#RUN apt-get update && install apt-get -y liblapack3 libboost1.62-dev gawk
RUN apt-get update && apt-get install -y pigz unzip libparallel-forkmanager-perl

# libparallel-forkmanager-perl required by iREAD 0.8.0
# argparse is requierd by iREAD 0.8.0

WORKDIR /opt

RUN conda create -n organoids python=2.7 &&\
    /bin/bash -c "source activate organoids && \
                  conda install -c bioconda deeptools=2.5.7 \
                                            rseqc=2.6.4 \
                                            samtools=1.8 \
                                            bedops=2.4.35 \
                                            argparse=1.4.0 \
                                            fastqc=0.11.7 \
                                            sambamba=0.6.6 \
                                            star=2.5.3a \
                                            subread=1.6.2 \
					    sra-tools=2.91 &&\
                 conda clean --all"

RUN curl -L -o BBMap_38.16.tar.gz "https://sourceforge.net/projects/bbmap/files/BBMap_38.16.tar.gz/download" &&\
    tar xzvf BBMap_38.16.tar.gz &&\
    rm BBMap_38.16.tar.gz

RUN curl -L -o iREAD_0.8.0.zip "http://www.genemine.org/codes/iREAD_0.8.0.zip" &&\
    (unzip iREAD_0.8.0.zip; exit 0) &&\
    ( find ./iREAD_0.8.0 -type f -name *.py -exec sed -i -e 's:/usr/bin/python:/usr/bin/env python:g' {} \; ) &&\
    ( find ./iREAD_0.8.0 -type f -name *.pl -exec sed -i -e 's:/usr/bin/perl:/usr/bin/env perl:g' {} \; ) &&\
    rm iREAD_0.8.0.zip

# This activate the organoids automatically, I guess there is a better way to do this
RUN echo "conda activate organoids" >> ~/.bashrc

ENV PATH=/opt/bbmap:/opt/iREAD_0.8.0:${PATH}
