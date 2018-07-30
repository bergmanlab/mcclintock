#!/bin/bash

# Establish variables
sample=$1
referencefolder=$2
te_cov_dir=$3
fastq1=$4
fastq2=$5
reference_genome=$6
consensus_te_seqs=$7
mcclintock_location=$8
processors=$9

# create tmp folder for intermediate files
tmp_dir=$referencefolder/te_cov_tmp
mkdir -p $tmp_dir
rm -rf $tmp_dir/*

# Hard mask reference genome to exclude reference TE sequences using repeatmasker
if [[ ! -f $reference_genome".masked" || ! -f $reference_genome".out.gff" ]]
then
    time RepeatMasker -dir $referencefolder -s -nolow -no_is -a $reference_genome -lib $consensus_te_seqs -pa $processors
    # RepeatMasker appears to override the custom database names during the ProcessRepeats step so this changes them back for
    # Drosophila, more rules like this may be needed for other reference genomes
    sed "s/McClintock-int/McClintock/g" $reference_genome".out.gff" > $referencefolder/tmp
    sed "s/POGON1/pogo/g" $referencefolder/tmp > $reference_genome".out.gff"
    perl $mcclintock_location/scripts/fixfastalinelength.pl $reference_genome".masked" 80 $reference_genome".masked2"
    mv $reference_genome".masked2" $reference_genome".masked"
fi

# create augmented genome by including TE sequences at the end of the reference genome.
ref_masked=$reference_genome".masked"
ref_augmented=$ref_masked".aug"
cat $ref_masked $consensus_te_seqs > $ref_augmented
samtools faidx $ref_augmented
bwa index $ref_augmented

# Map reads to the augmented genome.
bwa mem -t $processors -v 0 -R '@RG\tID:'$sample'\tSM:'$sample $ref_augmented $fastq1 $fastq2 > $tmp_dir/$sample.sam
sam=$tmp_dir/$sample.sam
samtools view -Sb -t $ref_augmented".fai" $sam | samtools sort - $tmp_dir/$sample
bam=$tmp_dir/$sample.bam
samtools index $bam

# get non-TE regions of the reference genome in bed format
rep_out=$reference_genome".out"
rep_bed=$reference_genome".bed"
rep_bed_sorted=$reference_genome".sorted.bed"
awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,".",strand}}' $rep_out > $rep_bed
bedtools sort -i $rep_bed > $rep_bed_sorted
grep chr $ref_augmented.fai | cut -f1,2 | sort > $referencefolder/$sample".sorted.genome"
bedtools slop -i $rep_bed_sorted -g $referencefolder/$sample".sorted.genome" -l 1 -r 0 > $referencefolder/$sample".fasta.out.zero.bed"
bedtools complement -i $referencefolder/$sample".fasta.out.zero.bed" -g $referencefolder/$sample".sorted.genome" | grep -v chrM > $referencefolder/$sample".fasta.out.complement.bed"
bed_nonte=$referencefolder/$sample".fasta.out.complement.bed"

# Calculate depth of coverage for every TE using Samtools depth module on the BAM file and normalize by average coverage in no TE regione of the normal reference genome.
grep ">" $consensus_te_seqs | sed 's/>//g' > $referencefolder/te_list
te_list=$referencefolder/te_list

genome_avg_depth=`samtools depth -b  $bed_nonte $bam | awk '{ total += $3 } END { print total/NR }'`
echo $genome_avg_depth > $referencefolder/genome_avg_depth
printf '%s\n' "TE Family" "Normalized Depth" | paste -sd ',' > $te_cov_dir/te_depth.csv
for te in `cat $te_list`
do
    te_depth=`samtools depth -r $te $bam | awk '{ total += $3 } END { print total/NR }'`
    te_depth_normalized=$(echo "$te_depth / $genome_avg_depth" | bc -l )
    printf '%s,%.2f\n' "$te" "$te_depth_normalized" | paste -sd ',' >> $te_cov_dir/te_depth.csv
done

# remove tmp folder
# rm -rf $tmp_dir