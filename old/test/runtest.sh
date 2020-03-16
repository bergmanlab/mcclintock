#!/bin/bash -l

run_dir=`pwd`

# Download fastq sequencing data from EBI
printf "\nDownloading and extracting fastq files...\n\n"

wget -nc -q ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_1.fastq.gz
wget -nc -q ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_2.fastq.gz

gunzip SRR800842_1.fastq.gz
gunzip SRR800842_2.fastq.gz

# Download the reference genome from UCSC (allows easy browsing of results)
printf "Downloading reference genome...\n\n"

wget -nc -q http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz
tar xvzf chromFa.tar.gz 
rm chromFa.tar.gz
# Combine the chromosomes together
cat chr*fa 2micron.fa > sacCer2.fasta
rm chr*fa 2micron.fa

# Download gff locations of reference TE copies
wget -nc -q http://files.figshare.com/287395/File_S2.txt
awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1;}' File_S2.txt > tmp
sed '1d;$d' tmp > reference_TE_locations.gff
rm File_S2.txt
rm tmp

# The TE families file and consensus TE fasta file are included in this folder

# Run the pipeline 
cd ..
bash mcclintock.sh -r $run_dir/sacCer2.fasta -c $run_dir/sac_cer_TE_seqs.fasta -g $run_dir/reference_TE_locations.gff -t $run_dir/sac_cer_te_families.tsv -1 $run_dir/SRR800842_1.fastq -2 $run_dir/SRR800842_2.fastq -p 4
