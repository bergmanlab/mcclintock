#!/usr/bin/bash

if (( $# > 0 ))
then
	mkdir $5
	# Run the TE locate pipeline.
	perl TE_locate.pl $4 $1 $3 $2 $5/ 600 3 1
	# Extract the relevant data from the output file 
	sed -e '1,2d' $5/"_600_reads3_acc1.info" > tmpfile
	
	#Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$5"_TE-locate"\" description=\"$5"_TE-locate"\"" > $5/$5_"telocate.bed"
	awk -F'[\t/]' '{printf $1"\t"$2"\t"$2"\t"$5; if($17=="old") printf "_"$17"\t0"; else if($17=="new") printf "_"$17"\t0"; if($14=="parallel") print "\t+"; else if($14 =="uncertain") print "\t."; else print "\t-";}' tmpfile >> $5/$5"_telocate_presort.bed"
	bedtools sort -i $5/$5_"telocate_presort.bed" >> $5/$5_"telocate.bed"
	rm tmpfile $5/$5"_telocate_presort.bed"
else
	echo "Supply lexically sorted SAM containing folder as option 1 (sort --temporary-directory=. <sam file> > <sorted sam file>)"
	echo "Supply fasta reference file as option 2"
	echo "Supply gff3 annotation of TEs inserts in reference as option 3"
	echo "Supply maximum memory (GB) for the software to use as option 4"
	echo "Supply output folder name as option 5"
fi

