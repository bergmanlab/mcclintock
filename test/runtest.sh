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
# Combine the chromosomes together with the TE sequences
cat chr*fa 2micron.fa sac_cer_TE_seqs.fa > sacCer2.fa
rm chr*fa 2micron.fa

# Download gff locations of reference TE copies
wget -nc -q http://files.figshare.com/287395/File_S2.txt
awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1;}' File_S2.txt > tmp
sed '1d;$d' tmp > reference_TE_locations.gff
rm File_S2.txt
rm tmp

# Add the locations of the sequences of the consensus TEs to the genome annotation
awk -F">" '/^>/ {if (seqlen){print seqlen}; printf $2"\t" ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' sac_cer_TE_seqs.fa > TE-lengths
while read TE length
do
    echo -e "$TE\treannotate\ttransposable_element\t1\t$length\t.\t+\t.\tID=$TE" >> reference_TE_locations.gff
done < TE-lengths

# The TE families file and consensus TE fasta file are included in this folder

# Run the pipeline 
cd ..
bash mcclintock.sh -r $run_dir/sacCer2.fa -c $run_dir/sac_cer_TE_seqs.fa -g $run_dir/reference_TE_locations.gff -t $run_dir/sac_cer_te_families.tsv -1 $run_dir/SRR800842_1.fastq -2 $run_dir/SRR800842_2.fastq -p 1
