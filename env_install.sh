#!/bin/bash

# remove env is already exist
conda remove -y -n MCCLINTOCK_main --all
conda remove -y -n MCCLINTOCK_cov --all
conda remove -y -n MCCLINTOCK_temp --all
conda remove -y -n MCCLINTOCK_ngstemapper --all
conda remove -y -n MCCLINTOCK_relocate --all
conda remove -y -n MCCLINTOCK_telocate --all
conda remove -y -n MCCLINTOCK_populationte --all
conda remove -y -n MCCLINTOCK_retroseq --all

# create specific conda envs for McC main pipeilne
conda create -y -n MCCLINTOCK_main
conda install -n MCCLINTOCK_main -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_main -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_main -y bwa=0.7.4
conda install -n MCCLINTOCK_main -y bedtools=2.17.0
conda install -n MCCLINTOCK_main -y bowtie=1.0.0
conda install -n MCCLINTOCK_main -y ucsc-blat=366
conda install -n MCCLINTOCK_main -y samtools=0.1.19
conda install -n MCCLINTOCK_main -y r=3.5.0
conda install -n MCCLINTOCK_main -y fastqc=0.11.2
conda install -n MCCLINTOCK_main -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_main -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_main -y bcftools=1.2
conda install -n MCCLINTOCK_main -y exonerate=2.4.0

# create specific conda envs for coverage module
conda create -y -n MCCLINTOCK_cov
conda install -n MCCLINTOCK_cov -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_cov -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_cov -y bwa=0.7.4
conda install -n MCCLINTOCK_cov -y bedtools=2.17.0
conda install -n MCCLINTOCK_cov -y bowtie=1.0.0
conda install -n MCCLINTOCK_cov -y ucsc-blat=366
conda install -n MCCLINTOCK_cov -y samtools=1.6
conda install -n MCCLINTOCK_cov -y r=3.5.0
conda install -n MCCLINTOCK_cov -y fastqc=0.11.2
conda install -n MCCLINTOCK_cov -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_cov -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_cov -y bcftools=1.2
conda install -n MCCLINTOCK_cov -y exonerate=2.4.0

# create specific conda envs for ngs_te_mapper
conda create -y -n MCCLINTOCK_ngstemapper
conda install -n MCCLINTOCK_ngstemapper -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_ngstemapper -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_ngstemapper -y bwa=0.7.4
conda install -n MCCLINTOCK_ngstemapper -y bedtools=2.17.0
conda install -n MCCLINTOCK_ngstemapper -y bowtie=1.0.0
conda install -n MCCLINTOCK_ngstemapper -y ucsc-blat=366
conda install -n MCCLINTOCK_ngstemapper -y samtools=0.1.19
conda install -n MCCLINTOCK_ngstemapper -y r=3.5.0
conda install -n MCCLINTOCK_ngstemapper -y fastqc=0.11.2
conda install -n MCCLINTOCK_ngstemapper -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_ngstemapper -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_ngstemapper -y bcftools=1.2
conda install -n MCCLINTOCK_ngstemapper -y exonerate=2.4.0

# create specific conda envs for temp
conda create -y -n MCCLINTOCK_temp
conda install -n MCCLINTOCK_temp -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_temp -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_temp -y bwa=0.7.4
conda install -n MCCLINTOCK_temp -y bedtools=2.17.0
conda install -n MCCLINTOCK_temp -y bowtie=1.0.0
conda install -n MCCLINTOCK_temp -y ucsc-blat=366
conda install -n MCCLINTOCK_temp -y samtools=0.1.19
conda install -n MCCLINTOCK_temp -y r=3.5.0
conda install -n MCCLINTOCK_temp -y fastqc=0.11.2
conda install -n MCCLINTOCK_temp -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_temp -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_temp -y bcftools=1.2
conda install -n MCCLINTOCK_temp -y exonerate=2.4.0

# create specific conda envs for relocate
conda create -y -n MCCLINTOCK_relocate
conda install -n MCCLINTOCK_relocate -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_relocate -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_relocate -y bwa=0.7.4
conda install -n MCCLINTOCK_relocate -y bedtools=2.17.0
conda install -n MCCLINTOCK_relocate -y bowtie=1.0.0
conda install -n MCCLINTOCK_relocate -y ucsc-blat=366
conda install -n MCCLINTOCK_relocate -y samtools=0.1.19
conda install -n MCCLINTOCK_relocate -y r=3.5.0
conda install -n MCCLINTOCK_relocate -y fastqc=0.11.2
conda install -n MCCLINTOCK_relocate -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_relocate -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_relocate -y bcftools=1.2
conda install -n MCCLINTOCK_relocate -y exonerate=2.4.0

# create specific conda envs for telocate
conda create -y -n MCCLINTOCK_telocate
conda install -n MCCLINTOCK_telocate -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_telocate -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_telocate -y bwa=0.7.4
conda install -n MCCLINTOCK_telocate -y bedtools=2.17.0
conda install -n MCCLINTOCK_telocate -y bowtie=1.0.0
conda install -n MCCLINTOCK_telocate -y ucsc-blat=366
conda install -n MCCLINTOCK_telocate -y samtools=0.1.19
conda install -n MCCLINTOCK_telocate -y r=3.5.0
conda install -n MCCLINTOCK_telocate -y fastqc=0.11.2
conda install -n MCCLINTOCK_telocate -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_telocate -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_telocate -y bcftools=1.2
conda install -n MCCLINTOCK_telocate -y exonerate=2.4.0

# create specific conda envs for MCCLINTOCK_populationte
conda create -y -n MCCLINTOCK_populationte
conda install -n MCCLINTOCK_populationte -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_populationte -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_populationte -y bwa=0.7.4
conda install -n MCCLINTOCK_populationte -y bedtools=2.17.0
conda install -n MCCLINTOCK_populationte -y bowtie=1.0.0
conda install -n MCCLINTOCK_populationte -y ucsc-blat=366
conda install -n MCCLINTOCK_populationte -y samtools=0.1.19
conda install -n MCCLINTOCK_populationte -y r=3.5.0
conda install -n MCCLINTOCK_populationte -y fastqc=0.11.2
conda install -n MCCLINTOCK_populationte -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_populationte -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_populationte -y bcftools=1.2
conda install -n MCCLINTOCK_populationte -y exonerate=2.4.0

# create specific conda envs for MCCLINTOCK_retroseq
conda create -y -n MCCLINTOCK_retroseq
conda install -n MCCLINTOCK_retroseq -y repeatmasker=4.0.7
conda install -n MCCLINTOCK_retroseq -y perl-bioperl-run=1.006900
conda install -n MCCLINTOCK_retroseq -y bwa=0.7.4
conda install -n MCCLINTOCK_retroseq -y bedtools=2.17.0
conda install -n MCCLINTOCK_retroseq -y bowtie=1.0.0
conda install -n MCCLINTOCK_retroseq -y ucsc-blat=366
conda install -n MCCLINTOCK_retroseq -y samtools=0.1.19
conda install -n MCCLINTOCK_retroseq -y r=3.5.0
conda install -n MCCLINTOCK_retroseq -y fastqc=0.11.2
conda install -n MCCLINTOCK_retroseq -y ucsc-fatotwobit=366
conda install -n MCCLINTOCK_retroseq -y ucsc-twobittofa=366
conda install -n MCCLINTOCK_retroseq -y bcftools=1.2
conda install -n MCCLINTOCK_retroseq -y exonerate=2.4.0