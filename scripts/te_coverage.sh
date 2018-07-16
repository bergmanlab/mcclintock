#!/bin/bash
# create tmp folder for intermediate files
tmp_dir=$referencefolder/te_cov_tmp
mkdir -p $tmp_dir
rm -rf $tmp_dir/*

# Hard mask reference genome to exclude reference TE sequences using repeatmasker
if [[ ! -f $reference_genome".masked" || ! -f $reference_genome".out.gff" ]]
then
    time RepeatMasker -dir $referencefolder -s -nolow -no_is -a $reference_genome -lib $consensus_te_seqs -pa $processors
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
grep chr $ref_augmented.fai | cut -f1,2 | sort > $referencefolder/dm6.sorted.genome
bedtools slop -i $rep_bed_sorted -g $referencefolder/dm6.sorted.genome -l 1 -r 0 > $referencefolder/dm6.fasta.out.zero.bed
bedtools complement -i $referencefolder/dm6.fasta.out.zero.bed -g $referencefolder/dm6.sorted.genome | grep -v chrM > $referencefolder/dm6.fasta.out.complement.bed
bed_nonte=$referencefolder/dm6.fasta.out.complement.bed

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