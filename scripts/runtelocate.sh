#!/usr/bin/bash

if (( $# > 0 ))
then

	# Establish variables
	sam_file=$1
	reference_genome=$2
	gff_te_locations=$3
	max_memory=$4
	sample=$5
	distance="$(($6 * 5))"
	outputfolder=$7

	mkdir -p $outputfolder/TE-locate

	# Run the TE locate pipeline.
	perl TE_locate.pl $max_memory $sam_file $gff_te_locations $reference_genome $outputfolder/TE-locate/ $distance 3 1
	# Extract the relevant data from the output file 
	sed -e '1,2d' $outputfolder/TE-locate/"_"$distance"_reads3_acc1.info" > $outputfolder/TE-locate/$sample"tmpfile"
	
	#Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$sample"_TE-locate"\" description=\"$sample"_TE-locate"\"" > $outputfolder/TE-locate/$sample"_telocate_redundant.bed"
	echo -e "track name=\"$sample"_TE-locate"\" description=\"$sample"_TE-locate"\"" > $outputfolder/TE-locate/$sample"_telocate_nonredundant.bed"

	awk -F'[\t/]' -v sample=$sample '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$8"\t"$5"_reference_"sample"_telocate_rp_\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$8"\t"$5"_non-reference_"sample"_telocate_rp_\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' $outputfolder/TE-locate/$sample"tmpfile" > $outputfolder/TE-locate/$sample"_telocate_presort.txt"
	bedtools sort -i $outputfolder/TE-locate/$sample"_telocate_presort.txt" > $outputfolder/TE-locate/$sample"_telocate_sorted_redundant.txt"

	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' $outputfolder/TE-locate/$sample"_telocate_sorted_redundant.txt" > $outputfolder/TE-locate/$sample"_telocate_counted_redundant.txt"
	cut -f1-3,5-7 $outputfolder/TE-locate/$sample"_telocate_counted_redundant.txt" >> $outputfolder/TE-locate/$sample"_telocate_redundant.bed"

	# Filter redundant predictions
	sort -k1,3 -k4rn $outputfolder/TE-locate/$sample"_telocate_counted_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $outputfolder/TE-locate/$sample"tmpfile"
	bedtools sort -i $outputfolder/TE-locate/$sample"tmpfile" >> $outputfolder/TE-locate/$sample"_telocate_nonredundant.bed"

	rm $outputfolder/TE-locate/$sample"tmpfile" $outputfolder/TE-locate/$sample"_telocate_presort.txt" $outputfolder/TE-locate/$sample"_telocate_counted_redundant.txt" $outputfolder/TE-locate/$sample"_telocate_sorted_redundant.txt"

else
	echo "Supply lexically sorted SAM containing folder as option 1 (sort --temporary-directory=. <sam file> > <sorted sam file>)"
	echo "Supply fasta reference file as option 2"
	echo "Supply gff3 annotation of TEs inserts in reference as option 3"
	echo "Supply maximum memory (GB) for the software to use as option 4"
	echo "Supply sample name as option 5"
	echo "Supply library insert size as option 6"
	echo "Supply output folder as option 7"
fi

