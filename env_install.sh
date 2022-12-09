#!/bin/bash
mamba remove -y -n MCCLINTOCK_main --all
mamba remove -y -n MCCLINTOCK_cov --all
mamba remove -y -n MCCLINTOCK_temp --all
mamba remove -y -n MCCLINTOCK_ngstemapper --all
mamba remove -y -n MCCLINTOCK_relocate --all
mamba remove -y -n MCCLINTOCK_telocate --all
mamba remove -y -n MCCLINTOCK_popoolationte --all
mamba remove -y -n MCCLINTOCK_retroseq --all

# create specific conda envs for McC main pipeilne
mamba create -y -n MCCLINTOCK_main python=3.6
mamba install -n MCCLINTOCK_main -y repeatmasker=4.0.7 perl-bioperl-run=1.006900 bwa=0.7.4 bedtools=2.17.0 samtools=0.1.19 fastqc=0.11.2 ucsc-fatotwobit=366 unzip=6.0 patch=2.7.6

# create specific conda envs for coverage module
mamba create -y -n MCCLINTOCK_cov python=3.6
mamba install -n MCCLINTOCK_cov -y repeatmasker=4.0.7 perl-bioperl-run=1.006900 bwa=0.7.4 bedtools=2.17.0 samtools=1.9 python=3.6 matplotlib

# create specific conda envs for ngs_te_mapper
mamba create -y -n MCCLINTOCK_ngstemapper python=3.6
mamba install -n MCCLINTOCK_ngstemapper -y bwa=0.7.4 bedtools=2.17.0 r=3.5.0

# create specific conda envs for temp
mamba create -y -n MCCLINTOCK_temp python=3.6
mamba install -n MCCLINTOCK_temp -y perl-bioperl-run=1.006900 bwa=0.7.4 bedtools=2.17.0 samtools=0.1.19 ucsc-fatotwobit=366 ucsc-twobittofa=366

# create specific conda envs for relocate
mamba create -y -n MCCLINTOCK_relocate python=3.6
mamba install -n MCCLINTOCK_relocate -y perl-bioperl-run=1.006900 bedtools=2.17.0 bowtie=1.0.0 ucsc-blat=366 samtools=0.1.19

# create specific conda envs for telocate
mamba create -y -n MCCLINTOCK_telocate python=3.6
mamba install -n MCCLINTOCK_telocate -y perl-bioperl-run=1.006900 bwa=0.7.4 bedtools=2.17.0 java-jdk

# create specific conda envs for MCCLINTOCK_populationte
mamba create -y -n MCCLINTOCK_popoolationte python=3.6
mamba install -n MCCLINTOCK_popoolationte -y perl-bioperl-run=1.006900 bwa=0.7.4 bedtools=2.17.0 samtools=0.1.19

# create specific conda envs for MCCLINTOCK_retroseq
mamba create -y -n MCCLINTOCK_retroseq python=3.6
mamba install -n MCCLINTOCK_retroseq -y perl-bioperl-run=1.006900 bwa=0.7.4 bedtools=2.17.0 samtools=0.1.19 bcftools=1.2 exonerate=2.4.0

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate MCCLINTOCK_cov
dir=`which openssl`
lib_dir=`echo $dir | sed 's/bin\/openssl/lib/'`
ln -s $lib_dir/libcrypto.so.1.1 $lib_dir/libcrypto.so.1.0.0
conda deactivate
