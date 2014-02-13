#!/bin/bash -l

test_dir=`pwd`

# Download fastq sequencing data from EBI

printf "\nDownloading and extracting fastq files...\n\n"

wget -nc ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_1.fastq.gz
wget -nc ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_2.fastq.gz

gunzip SRR800842_1.fastq.gz
gunzip SRR800842_2.fastq.gz

# Download the reference genome from UCSC (allows easy browsing of results)

printf "Downloading reference genome...\n\n"

wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz
tar xvzf chromFa.tar.gz 
rm chromFa.tar.gz
cat chr*fa 2micron.fa > sacCer2.fa
rm chr*fa 2micron.fa

# Download gff locations of reference TE copies
wget -nc http://files.figshare.com/287395/File_S2.txt
awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1;}' File_S2.txt > tmp
sed '1d;$d' tmp > reference_TE_locations.gff
rm File_S2.txt
rm tmp

# The TE families file and consensus TE fasta file are included in this folder

# Run the pipeline 
cd ..
bash mcclintock $test_dir/sacCer2.fa $test_dir/sac_cer_TE_seqs.fa $test_dir/reference_TE_locations.gff $test_dir/te_families $test_dir/SRR800842_1.fastq $test_dir/SRR800842_2.fastq
