#!/bin/bash

if (( $# > 0 ))
then
	reference=${3##*/}
	reference=${reference%%.*}
	perl splitforRetroSeq.pl $1 $reference
	samplename=`basename $2 .bam`
	mkdir $samplename
	bin/retroseq.pl -discover -bam $2 -eref $reference"elementlist" -output $samplename/$samplename.discovery 
	bin/retroseq.pl -call -bam $2 -input $samplename/$samplename.discovery -ref $3 -output $samplename/$samplename.calling -orientate yes
	
	# Extract the relevant results
	echo -e "track name=\"$samplename"_RetroSeq"\" description=\"$samplename"_RetroSeq"\"" > $samplename/$samplename"_retroseq.bed"
	awk '$1!~/#/{print $0}' $samplename/$samplename.calling.PE.vcf >> tmp
	awk -F'[=,\t:]' '{if ($20 >= 6) print $1"\t"$11"\t"$12"\t"$10"_new\t0\t."}' tmp >> $samplename/$samplename"_retroseq_presort.bed"
	bedtools sort -i $samplename/$samplename"_retroseq_presort.bed" >> $samplename/$samplename"_retroseq.bed"
	rm tmp $samplename/$samplename"_retroseq_presort.bed"
	
else
	echo "Supply TE database as option 1"
	echo "Supply BAM file as option 2"
	echo "Supply fasta reference genome used to make BAM as option 3"
fi
