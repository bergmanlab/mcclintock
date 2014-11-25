#!/usr/bin/bash

if (( $# > 0 ))
then
    sample=$5
	mkdir $sample
	distance="$(($6 * 5))"
	# Run the TE locate pipeline.
	perl TE_locate.pl $4 $1 $3 $2 $sample/ $distance 3 1
	# Extract the relevant data from the output file 
	sed -e '1,2d' $sample/"_"$distance"_reads3_acc1.info" > $sample/$sample"tmpfile"
	
	#Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$sample"_TE-locate"\" description=\"$sample"_TE-locate"\"" > $sample/$sample"_telocate_redundant.bed"
	echo -e "track name=\"$sample"_TE-locate"\" description=\"$sample"_TE-locate"\"" > $sample/$sample"_telocate_nonredundant.bed"



	awk -F'[\t/]' -v sample=$sample '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$8"\t"$5"_"$17"_"sample"_telocate_rp_\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$8"\t"$5"_"$17"_"sample"_telocate_rp_\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' $sample/$sample"tmpfile" > $sample/$sample"_telocate_presort.txt"
	bedtools sort -i $sample/$sample"_telocate_presort.txt" > $sample/$sample"_telocate_sorted_redundant.txt"

	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' $sample/$sample"_telocate_sorted_redundant.txt" > $sample/$sample"_telocate_counted_redundant.txt"
	cut -f1-3,5-7 $sample/$sample"_telocate_counted_redundant.txt" >> $sample/$sample"_telocate_redundant.bed"

	# Filter redundant predictions
	sort -k1,3 -k4rn $sample/$sample"_telocate_counted_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $sample/$sample"tmpfile"
	bedtools sort -i $sample/$sample"tmpfile" >> $sample/$sample"_telocate_nonredundant.bed"

	rm $sample/$sample"tmpfile" $sample/$sample"_telocate_presort.txt" $sample/$sample"_telocate_counted_redundant.txt" $sample/$sample"_telocate_sorted_redundant.txt"

else
	echo "Supply lexically sorted SAM containing folder as option 1 (sort --temporary-directory=. <sam file> > <sorted sam file>)"
	echo "Supply fasta reference file as option 2"
	echo "Supply gff3 annotation of TEs inserts in reference as option 3"
	echo "Supply maximum memory (GB) for the software to use as option 4"
	echo "Supply output folder name as option 5"
	echo "Supply library insert size as option 6"
fi

