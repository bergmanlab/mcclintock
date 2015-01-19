#!/bin/bash

if (( $# > 0 ))
then

	# Establish variables
	test_dir=`pwd`
	bam=$1
	sam=$2
	consensus_te_seqs=$3
	bed_te_locations_file=$4
	te_families=$5
	median_insertsize=$6
	sample=$7
	processors=$8

	mkdir $sample
	chmod 755 $sample

	if [[ $median_insertsize == "none" ]]
	then
		median_insertsize=`cut -f9 $sam | sort -n | awk '{if ($1 > 0) ins[reads++]=$1; } END { print ins[int(reads/2)]; }'`
	fi

	# Create soft link to differently named files for TEMP
	cp -s $bam $sample/$sample.sorted.bam
	cp -s $bam.bai $sample/$sample.sorted.bam.bai

	bash $test_dir/scripts/TEMP_Insertion.sh -x 30 -i $test_dir/$sample/$sample.sorted.bam -s $test_dir/scripts -r $consensus_te_seqs -t $bed_te_locations_file -u $te_families -m 1 -f $median_insertsize -c $processors -o $test_dir/$sample

	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $sample/$sample"_temp_raw.bed"
	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $sample/$sample"_temp_redundant.bed"
	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $sample/$sample"_temp_nonredundant.bed"

	sed '1d' $sample/$sample".insertion.refined.bp.summary" | awk '{if ($4 == "sense" || $4 == "antisense"); else print $0}' | awk -v sample=$sample '{ printf $1"\t"$9"\t"$11"\t"$7"\t"$4"_new_"sample"_temp_\t0\t"; if ( $5 == "sense" ) printf "+"; else printf "-"; print "\t"$6"\t"$10"\t"$12}' > $sample/$sample"_temp_presort_raw.txt"

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
	echo "Supply a bam file created by bwa mem as option 1"
	echo "Supply a sam file created by bwa mem as option 2"
	echo "Supply a consensus TE sequence file (fasta) as option 3"
	echo "Supply reference insertion locations in bed format as option 4"
	echo "Supply TE hierarchy file as option 5"
	echo "Supply the median insert length as option 6 or 'none' if it is not known"
	echo "Supply the sample name as option 7"
	echo "Supply the number of processors to use as option 8"
fi
