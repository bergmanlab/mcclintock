#!/bin/bash

if (( $# == 10 ))
then

	# Establish variables
	test_dir=`pwd`
	bam=$1
	reference2bit=$2
	consensus_te_seqs=$3
	bed_te_locations_file=$4
	te_families=$5
	median_insertsize=$6
	sample=$7
	processors=$8
	outputfolder=$9
	gff_te_locations=${10}

	mkdir $outputfolder/TEMP
	chmod 755 $outputfolder/ $outputfolder/TEMP

	# Create soft link to differently named files for TEMP
	cp -s $bam $outputfolder/TEMP/$sample.sorted.bam
	cp -s $bam.bai $outputfolder/TEMP/$sample.sorted.bam.bai

	# Run presence module
	bash $test_dir/scripts/TEMP_Insertion.sh -x 30 -i $outputfolder/TEMP/$sample.sorted.bam -s $test_dir/scripts -r $consensus_te_seqs -t $bed_te_locations_file -u $te_families -m 1 -f $median_insertsize -c $processors -o $outputfolder/TEMP

	# Run absence module - currently saved in originalmethodresults for McClintock but not parsed anywhere
	bash $test_dir/scripts/TEMP_Absence.sh -i $outputfolder/TEMP/$sample.sorted.bam -s $test_dir/scripts -r $bed_te_locations_file -t $reference2bit -f $median_insertsize -c $processors -o $outputfolder/TEMP

	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $outputfolder/TEMP/$sample"_temp_raw.bed"
	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $outputfolder/TEMP/$sample"_temp_redundant.bed"
	echo -e "track name=\"$sample"_TEMP"\" description=\"$sample"_TEMP"\"" > $outputfolder/TEMP/$sample"_temp_nonredundant.bed"

	# Deal with malformed predictions and predictions that have no TE call and begin to convert data
	awk '{if ($3>=$2 && $3 > 0 && $2 > 0) print $0}'  $outputfolder/TEMP/$sample".insertion.refined.bp.summary" | awk '{if ($4 == "sense" || $4 == "antisense"); else print $0}' | awk -v sample=$sample '{ printf $1"\t"$2"\t"$3"\t"$7"\t"$4"_non-reference_"$8"_"sample"_temp_\t0\t"; if ( $5 == "sense" ) printf "+"; else printf "-"; print "\t"$6"\t"$10"\t"$12"\t"$9"\t"$11"\t"$8}' > $outputfolder/TEMP/$sample"_temp_presort_raw.txt"

	# Take the absent TEs and inverse them to report those TEs with no evidence of absence
	cut -f1-3 $outputfolder/TEMP/$sample".absence.refined.bp.summary" | sed '1d' > $outputfolder/TEMP/$sample".absent.bed"

	# Add a false record in case of a situation where no absent TEs are detected
	echo -e 'empty\t0\t1' >> $outputfolder/TEMP/$sample".absent.bed"
	bedtools subtract -A -a $gff_te_locations -b $outputfolder/TEMP/$sample".absent.bed" | awk -F'[\t=;]' -v sample=$sample '{print $1"\t"$4-1"\t"$5"\t!\t"$10"_reference_"sample"_temp_nonab_\t0\t"$7"\t!\t!\t!\t!\t!\t!"}' >> $outputfolder/TEMP/$sample"_temp_presort_raw.txt"

	bedtools sort -i $outputfolder/TEMP/$sample"_temp_presort_raw.txt" > $outputfolder/TEMP/$sample"_temp_sorted_raw.txt"

	awk '{ if ($5~/_nonab_/) {print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7} else {if ($9 > 0 && $10 > 0) printf $1"\t"$11-1"\t"$12"\t"$4"\t"$5"sr_"; else printf $1"\t"$2-1"\t"$3"\t"$4"\t"$5"rp_"; print NR"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$13}}' $outputfolder/TEMP/$sample"_temp_sorted_raw.txt" > $outputfolder/TEMP/$sample"_temp_sorted_precut_raw.txt"

	cut -f1-3,5-7 $outputfolder/TEMP/$sample"_temp_sorted_precut_raw.txt" >> $outputfolder/TEMP/$sample"_temp_raw.bed"

	# Filter for insertions with support for both ends or not absent reference sequences
	awk '{if ($8 == "1p1" && $11 > 0.1) print $0; else if ($4 == "!") print $0}' $outputfolder/TEMP/$sample"_temp_sorted_precut_raw.txt" > $outputfolder/TEMP/$sample"_temp_sorted_redundant.txt"

	cut -f1-3,5-7 $outputfolder/TEMP/$sample"_temp_sorted_redundant.txt" >> $outputfolder/TEMP/$sample"_temp_redundant.bed"

	# Filter out redundant insertions
	sort -k1,3 -k4rn $outputfolder/TEMP/$sample"_temp_sorted_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $outputfolder/TEMP/tmp
	bedtools sort -i $outputfolder/TEMP/tmp >> $outputfolder/TEMP/$sample"_temp_nonredundant.bed"

	rm $outputfolder/TEMP/tmp $outputfolder/TEMP/$sample"_temp_sorted_redundant.txt" $outputfolder/TEMP/$sample"_temp_presort_raw.txt" $outputfolder/TEMP/$sample"_temp_sorted_raw.txt" $outputfolder/TEMP/$sample"_temp_sorted_precut_raw.txt"

else
	echo "Supply a bam file created by bwa mem as option 1"
	echo "Supply a reference 2bit file as option 2"
	echo "Supply a consensus TE sequence file (fasta) as option 3"
	echo "Supply reference insertion locations in bed format as option 4"
	echo "Supply TE hierarchy file as option 5"
	echo "Supply the median insert length as option 6"
	echo "Supply the sample name as option 7"
	echo "Supply the number of processors to use as option 8"
	echo "Supply output folder as option 9"
	echo "Supply reference insertion locations in gff format with correct family names as the ID as option 10"
fi
