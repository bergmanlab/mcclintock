#!/bin/bash

cd $(dirname "$0")

if [[ ! -f sacCer2.fasta ]]
then
    printf "Downloading reference genome...\n\n"
    wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz
    tar xvzf chromFa.tar.gz
    rm chromFa.tar.gz
    cat chr*fa 2micron.fa > sacCer2.fasta
    rm chr*fa 2micron.fa
fi

if [[ ! -f SRR800842_1.fastq ]]
then
    printf "\nDownloading and extracting fastq file 1...\n\n"
    wget -nc ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_1.fastq.gz
    gunzip SRR800842_1.fastq.gz
fi

if [[ ! -f SRR800842_2.fastq ]]
then
    printf "\nDownloading and extracting fastq file 2...\n\n"
    wget -nc ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_2.fastq.gz
    gunzip SRR800842_2.fastq.gz
fi