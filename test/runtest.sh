#!/bin/bash -l

# This script will download all the data necessary to run all methods in 
# the pipeline and will launch the individual methods. It can be edited
# for use with different samples or organisms. (AS LONG AS CERTAIN INFORMATION IS KNOWN? TSDS?)

# Load modules necessary to run the software on kadmon
source /etc/profile.d/modules.sh
module load apps/bwa/0.7.5a
module load apps/perl/5.18.1
module load apps/RepeatMasker/3.0.3-p1
module load apps/samtools/0.1.19

testdir=`pwd`

# Set up folder structure
mkdir $testdir/fastq
mkdir $testdir/bam
mkdir $testdir/sam

# Download fastq sequencing data from EBI
wget -nc ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_1.fastq.gz
wget -nc ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_2.fastq.gz

printf "Extracting fastq files...\n\n"

gunzip SRR800842_1.fastq.gz
gunzip SRR800842_2.fastq.gz
mv SRR800842_1.fastq $testdir/fastq
fastq1=$testdir/fastq/SRR800842_1.fastq
mv SRR800842_2.fastq $testdir/fastq
fastq2=$testdir/fastq/SRR800842_2.fastq
samplename=SRR800842

# Download the reference genome from UCSC (allows easy browsing of results)

printf "Downloading reference genome...\n\n"

wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz
tar xvzf chromFa.tar.gz 
rm chromFa.tar.gz
cat chr*fa 2micron.fa > sacCer2.fa
rm chr*fa 2micron.fa
mv sacCer2.fa reference
reference_genome=$testdir/reference/sacCer2.fa
samtools faidx $reference_genome
bwa index $reference_genome

# Download gff locations of reference TE copies
wget -nc http://files.figshare.com/287395/File_S2.txt
awk '{print $3"\treannotate\ttransposable_element\t"$4"\t"$5"\t.\t"$6"\t.\tID="$1";Name="$1";Alias="$1}' File_S2.txt > tmp
sed '1d;$d' tmp > $testdir/reference/reference_TE_locations.gff
rm File_S2.txt
rm tmp
reference_tes=$testdir/reference/reference_TE_locations.gff

# Extract sequence of all reference TE copies
awk -F'[=\t]' 'BEGIN {OFS = "\t"}; {print $1,$2,$12,$4,$5,$6,$7,$8,$9"="$10"="$11"="$12}' $reference_tes > edited.gff
bedtools getfasta -name -fi $reference_genome -bed edited.gff -fo $testdir/reference/all_sac_Cer_TE_seqs.fa
rm edited.gff
all_TE_seqs=$testdir/reference/all_sac_Cer_TE_seqs.fa

# Set path for TE database
te_database=$testdir/reference/sac_cer_TE_seqs.fa

printf "Create bam alignment...\n\n"

# Create sam and bam files for input
bwa mem -v 1 $reference_genome $fastq1 $fastq2 > $testdir/bam/SRR800842.sam
sam=$testdir/bam/SRR800842.sam
sort --temporary-directory=. $sam > $testdir/bam/sortedSRR800842.sam
rm $sam
mv $testdir/bam/sortedSRR800842.sam $testdir/sam/SRR800842.sam 
sam=$testdir/sam/SRR800842.sam
sam_folder=$testdir/sam

samtools view -Sb $sam > $testdir/bam/SRR800842.bam
samtools sort $testdir/bam/SRR800842.bam $testdir/bam/sortedSRR800842
rm $testdir/bam/SRR800842.bam
mv $testdir/bam/sortedSRR800842.bam $testdir/bam/SRR800842.bam
bam=$testdir/bam/SRR800842.bam
samtools index $bam

sample_name=`basename $sam .sam`

# Run RelocaTE

printf "\nRunning RelocaTE pipeline...\n\n"

# Add TSD lengths to consensus TE sequences
awk '{if (/>/) print $0" TSD=UNK"; else print $0}' $te_database > $testdir/reference/relocaTE_sac_cer_TE_seqs.fa
relocate_te_database=$testdir/reference/relocaTE_sac_cer_TE_seqs.fa

cd ../RelocaTE
bash runrelocate.sh $relocate_te_database $reference_genome $testdir/fastq $sample_name $reference_tes

# Run ngs_te_mapper pipeline

printf "\nRunning ngs_te_mapper pipeline...\n\n"

fasta1=$testdir/fastq/SRR800842_1.fa
fasta2=$testdir/fastq/SRR800842_2.fa

cd ../ngs_te_mapper/

bash runngstemapper.sh $te_database $reference_genome $sample_name $fasta1 $fasta2 

# Run RetroSeq

printf "\nRunning RetroSeq pipeline...\n\n"

cd ../RetroSeq
bash runretroseq.sh $te_database $bam $reference_genome 

# Run TE-locate

printf "\nRunning TE-locate pipeline...\n\n"

# Create te_hierachy
printf "insert\tid\tfamily\tsuperfamily\tsuborder\torder\tclass\tproblem\n" > $testdir/reference/te_hierarchy
grep \> $all_TE_seqs | awk -F'[>_]' '{if ($0 ~ /TY3_1p/) {print $2"_"$3"_"$4"_"$5"\t"$4"_"$5"\t"$4"_"$5"\tGypsy\t LTR\tLTR\tRNA\t0"} else {{printf $2"_"$3"_"$4"\t"$4"\t"$4"\t"} if ($4 ~ /TY3/) {printf "Gypsy"} else {printf "Copia"} printf "\t LTR\tLTR\tRNA\t0\n"} }' >> $testdir/reference/te_hierarchy
te_hierarchy=$testdir/reference/te_hierarchy
# Adjust hierachy levels
cd ../TE-locate
sed '1d' $te_hierarchy | cut -f 1-2 > $testdir/reference/TE-locate_conversion
conversion_file=$testdir/reference/TE-locate_conversion
perl TE_hierarchy.pl $reference_tes $conversion_file Alias
telocate_reference_tes=${reference_tes%.*}
telocate_reference_tes=$telocate_reference_tes"_HL.gff"

bash runtelocate.sh $sam_folder $reference_genome $telocate_reference_tes 2 $sample_name

# Run PoPoolationTE

printf "\nRunning PoPoolationTE pipeline...\n\n"

cd ../popoolationte
bash runpopoolationte.sh $reference_genome $all_TE_seqs $te_hierarchy $fastq1 $fastq2 $reference_tes

