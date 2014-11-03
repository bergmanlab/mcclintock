#!/usr/bin/bash

if (( $# > 0 ))
then
	mkdir $5
	distance="$(($6 * 5))"
	# Run the TE locate pipeline.
	perl TE_locate.pl $4 $1 $3 $2 $5/ $distance 3 1
	# Extract the relevant data from the output file 
	sed -e '1,2d' $5/"_"$distance"_reads3_acc1.info" > $5/$5"tmpfile"
	
	#Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$5"_TE-locate"\" description=\"$5"_TE-locate"\"" > $5/$5"_telocate_duplicated.bed"
	echo -e "track name=\"$5"_TE-locate"\" description=\"$5"_TE-locate"\"" > $5/$5"_telocate.bed"
	awk -F'[\t/]' sample=$5 '{printf $1"\t"; if($17=="old") printf $2-1"\t"$2"\t"$5"_"$17"_"sample"_telocate_rp\t0"; else if($17=="new") printf $2-1"\t"$2"\t"$5"_"$17"_"sample"_telocate_rp\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' $5/$5"tmpfile" >> $5/$5"_telocate_presort.bed"

    sort -k1,3 -k4rn $5/$5"_telocate_presort.bed" | sort -u -k1,3 | cut -f1-3,5- > $5/$5"tmpfile"
    bedtools sort -i $5/$5"tmpfile" >> $5/$5"_telocate.bed"
	bedtools sort -i $5/$5"_telocate_presort.bed" >> $5/$5"_telocate_duplicated.bed"
	rm $5/$5"tmpfile" $5/$5"_telocate_presort.bed"
else
	echo "Supply lexically sorted SAM containing folder as option 1 (sort --temporary-directory=. <sam file> > <sorted sam file>)"
	echo "Supply fasta reference file as option 2"
	echo "Supply gff3 annotation of TEs inserts in reference as option 3"
	echo "Supply maximum memory (GB) for the software to use as option 4"
	echo "Supply output folder name as option 5"
	echo "Supply library insert size as option 6"
fi

