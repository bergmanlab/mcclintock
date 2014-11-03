#!/bin/bash

if (( $# > 0 ))
then
	reference=${3##*/}
	reference=${reference%%.*}

	if [ ! -d $reference ]
	then
		mkdir $reference
		# Create the individual input files required
		perl splitforRetroSeq.pl $1 $reference
		awk -F'[\t]' '{print $1}' $reference"_elementlist" > list
		awk -F".fa" '{print $1".bed"}' $reference"_elementlist" > $reference"_locationlist"

		# Link individual inserts to the bed file they must be annotated in
		# env LC_COLLATE=C is added as there was inconsistency between the sorting order when run on certain machines.
		# Errors from RetroSeq not finding TE specific bed files may originate here if these sorts fail.
		join -1 2 -2 1 -o 1.1,2.2 -t $'\t' <(env LC_COLLATE=C sort -k2 $5) <(env LC_COLLATE=C sort $reference"_locationlist") > tmp

		while read element file
		do
			grep -w $element $4 >> $file
		done < tmp

		rm tmp
		rm list
	fi

	samplename=`basename $2 .bam`
	mkdir $samplename

	# To run RetroSeq with the more computationally expensive Exonerate step uncomment the next line
	# and comment the alternate discovery step below it
	# bin/retroseq.pl -discover -bam $2 -eref $reference"_elementlist" -output $samplename/$samplename.discovery
	bin/retroseq.pl -discover -bam $2 -refTEs $reference"_locationlist" -output $samplename/$samplename.discovery

	bin/retroseq.pl -call -bam $2 -input $samplename/$samplename.discovery -filter $reference"_locationlist" -ref $3 -output $samplename/$samplename.calling -orientate yes

	# Extract the relevant results
	echo -e "track name=\"$samplename"_RetroSeq"\" description=\"$samplename"_RetroSeq"\"" > $samplename/$samplename"_retroseq_unfiltered.bed"
    echo -e "track name=\"$samplename"_RetroSeq"\" description=\"$samplename"_RetroSeq"\"" > $samplename/$samplename"_retroseq_duplicated.bed"
    echo -e "track name=\"$samplename"_RetroSeq"\" description=\"$samplename"_RetroSeq"\"" > $samplename/$samplename"_retroseq.bed"

	awk '$1!~/#/{print $0}' $samplename/$samplename.calling.PE.vcf >> $samplename/tmp
	awk -F'[=,\t:]' sample=$samplename '{if ($21 >= 6) print $1"\t"$11"\t"$12"\t"$10"_new_"sample"_retroseq_rp\t0\t."}' $samplename/tmp >> $samplename/$samplename"_retroseq_presort.bed"
    sort -k1,3 -k4rn $samplename/$samplename"_retroseq_presort.bed" | sort -u -k1,3 | cut -f1-3,5- > $samplename/tmp

	bedtools sort -i $samplename/$samplename"_retroseq_presort.bed" >> $samplename/$samplename"_retroseq_duplicated.bed"
    bedtools sort -i $samplename/tmp >> $samplename/$samplename"_retroseq.bed"
	rm $samplename/tmp $samplename/$samplename"_retroseq_presort.bed"

else
	echo "Supply TE database as option 1"
	echo "Supply BAM file as option 2"
	echo "Supply fasta reference genome used to make BAM as option 3"
	echo "Supply bed of TE family locations as option 4"
	echo "Supply TE hierarchy file as option 5"
fi
