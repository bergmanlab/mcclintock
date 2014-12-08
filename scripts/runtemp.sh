#!/bin/bash

if (( $# > 0 ))
then

    # Establish variables
	test_dir=`pwd`
	reference_genome=$1
	fastq1=$2
	fastq2=$3
	consensus_te_seqs=$4
	bed_te_locations_file=$5
	te_families=$6
	median_insertsize=$7
	sample=$8
	processors=$9

	mkdir $sample
	chmod 755 $sample

	# Have to use bwa aln for TEMP - this can be replaced with the commented lines below if fixed:
	bwa aln -t $processors $reference_genome $fastq1 > $sample/1.sai
	bwa aln -t $processors $reference_genome $fastq2 > $sample/2.sai

	bwa sampe $reference_genome $sample/1.sai $sample/2.sai $fastq1 $fastq2 | ../scripts/samtest.pl > $test_dir/$sample/$sample.sam

	if [[ $median_insertsize == "none" ]]
	then
		median_insertsize=`cut -f9 $sam | sort -n | awk '{if ($1 > 0) ins[reads++]=$1; } END { print ins[int(reads/2)]; }'`
	fi

	samtools view -Sb $test_dir/$sample/$sample.sam > $test_dir/$sample/$sample.bam
	rm $test_dir/$sample/$sample.sam

	# Get stats of bam file from samtools
	samtools flagstat $test_dir/$sample/$sample.bam > ../$genome/$sample/results/qualitycontrol/bwaaln_bamstats.txt

	samtools sort $test_dir/$sample/$sample.bam $test_dir/$sample/$sample.sorted
	rm $test_dir/$sample/$sample.bam
	samtools index $test_dir/$sample/$sample.sorted.bam

	# Lines to uncomment if TEMP switches to bwa mem - also need to move bam deletion to after TEMP, also resolve median insert size calc
	# Create soft link to differently named files for TEMP
	# cp -s $bam $sample/$sample.sorted.bam
	# cp -s $bam.bai $sample/$sample.sorted.bam.bai

	bash $test_dir/scripts/TEMP_Insertion.sh -i $test_dir/$sample/$sample.sorted.bam -s $test_dir/scripts -r $consensus_te_seqs -t $bed_te_locations_file -u $te_families -m 1 -f $median_insertsize -c $processors -o $test_dir/$sample

	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $sample/$sample"_temp_raw.bed"
	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $sample/$sample"_temp_redundant.bed"
	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $sample/$sample"_temp_nonredundant.bed"

	sed '1d' $sample/$sample.insertion.refined.bp.summary | awk '{if ($4 == "sense" || $4 == "antisense"); else print $0}' | awk -v sample=$sample '{ printf $1"\t"$9"\t"$11"\t"$7"\t"$4"_new_"sample"_temp_\t0\t"; if ( $5 == "sense" ) printf "+"; else printf "-"; print "\t"$6"\t"$10"\t"$12}' > $sample/$sample"_temp_presort_raw.txt"

	bedtools sort -i $sample/$sample"_temp_presort_raw.txt" > $sample/$sample"_temp_sorted_raw.txt"

	awk '{ printf $1"\t"$2"\t"$3"\t"$4"\t"$5; if ($9 > 0 && $10 > 0) printf "sr_"; else printf "rp_"; print NR"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $sample/$sample"_temp_sorted_raw.txt" > $sample/$sample"_temp_sorted_precut_raw.txt"

	cut -f1-3,5-7 $sample/$sample"_temp_sorted_precut_raw.txt" >> $sample/$sample"_temp_raw.bed"

	# Filter for insertions with support for both ends
	awk '{if ($8 == "1p1") print $0}' $sample/$sample"_temp_sorted_precut_raw.txt" > $sample/$sample"_temp_sorted_redundant.txt"

	cut -f1-3,5-7 $sample/$sample"_temp_sorted_redundant.txt" >> $sample/$sample"_temp_redundant.bed"

	# Filter out redundant insertions
	sort -k1,3 -k4rn $sample/$sample"_temp_sorted_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $sample/tmp
	bedtools sort -i $sample/tmp >> $sample/$sample"_temp_nonredundant.bed"

	rm $sample/tmp $sample/$sample"_temp_sorted_redundant.txt" $sample/$sample"_temp_presort_raw.txt" $sample/$sample"_temp_sorted_raw.txt" $sample/$sample"_temp_sorted_precut_raw.txt"

else
	echo "Supply fasta reference file as option 1"
	echo "Supply sequence fastas as options 2 and 3"
	echo "Supply a consensus TE sequence file (fasta) as option 4"
	echo "Supply reference insertion locations in bed format as option 5"
	echo "Supply TE hierarchy file as option 6"
	echo "Supply the median insert length as option 7"
	echo "Supply the sample name as option 8"
	echo "Supply the number of processors to use as option 9"
fi
