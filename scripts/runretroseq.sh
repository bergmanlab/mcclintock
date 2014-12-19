#!/bin/bash

if (( $# > 0 ))
then

    # Establish variables
    consensus_te_sequences=$1
    bam_file=$2
    reference_genome=$3
	reference=${3##*/}
	reference=${reference%%.*}
    bed_te_location=$4
    te_hierarchy=$5

	if [ ! -d $reference ]
	then
		mkdir $reference
		# Create the individual input files required
		perl splitforRetroSeq.pl $consensus_te_sequences $reference
		awk -F".fa" '{print $1".bed"}' $reference"_elementlist" > $reference"_locationlist"

		while read element file
		do
			awk -v element=$element '{if ($2==element) print $1}' $te_hierarchy > intfile
			while read insert
			do
				grep -w $insert $bed_te_location >> $file
			done < intfile
		done < $reference"_locationlist"
		rm intfile
	fi

	samplename=`basename $bam_file .bam`
	mkdir $samplename

	# To run RetroSeq with the more computationally expensive Exonerate step uncomment the next line
	# and comment the alternate discovery step below it
	# bin/retroseq.pl -discover -bam $bam_file -eref $reference"_elementlist" -output $samplename/$samplename.discovery
	bin/retroseq.pl -discover -bam $bam_file -refTEs $reference"_locationlist" -output $samplename/$samplename.discovery

	bin/retroseq.pl -call -bam $bam_file -input $samplename/$samplename.discovery -filter $reference"_locationlist" -ref $reference_genome -output $samplename/$samplename.calling -orientate yes

	# Extract the relevant results
	echo -e "track name=\"$samplename"_RetroSeq"\" description=\"$samplename"_RetroSeq"\"" > $samplename/$samplename"_retroseq_raw.bed"
	echo -e "track name=\"$samplename"_RetroSeq"\" description=\"$samplename"_RetroSeq"\"" > $samplename/$samplename"_retroseq_redundant.bed"
	echo -e "track name=\"$samplename"_RetroSeq"\" description=\"$samplename"_RetroSeq"\"" > $samplename/$samplename"_retroseq_nonredundant.bed"

	awk '$1!~/#/{print $0}' $samplename/$samplename.calling.PE.vcf > $samplename/tmp
	awk -F'[=,\t:]' -v sample=$samplename '{print $1"\t"$11"\t"$12"\t"$6"\t"$10"_new_"sample"_retroseq_rp_\t0\t.\t"$21}' $samplename/tmp > $samplename/$samplename"_retroseq_presort.txt"
	bedtools sort -i $samplename/$samplename"_retroseq_presort.txt" > $samplename/$samplename"_retroseq_sorted_raw.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8}' $samplename/$samplename"_retroseq_sorted_raw.txt" > $samplename/$samplename"_retroseq_sorted_counted_raw.txt"
	cut -f1-3,5-7 $samplename/$samplename"_retroseq_sorted_counted_raw.txt" >> $samplename/$samplename"_retroseq_raw.bed"

	# Filter results
	awk '{if ($8 >= 6) print $0}' $samplename/$samplename"_retroseq_sorted_counted_raw.txt" > $samplename/$samplename"_retroseq_redundant.txt"
	cut -f1-3,5-7 $samplename/$samplename"_retroseq_redundant.txt" >> $samplename/$samplename"_retroseq_redundant.bed"

	# Filter for redundant predictions
	sort -k1,3 -k4rn $samplename/$samplename"_retroseq_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $samplename/tmp
	bedtools sort -i $samplename/tmp >> $samplename/$samplename"_retroseq_nonredundant.bed"

	rm $samplename/tmp $samplename/$samplename"_retroseq_redundant.txt" $samplename/$samplename"_retroseq_sorted_counted_raw.txt" $samplename/$samplename"_retroseq_sorted_raw.txt" $samplename/$samplename"_retroseq_presort.txt"

else
	echo "Supply TE database as option 1"
	echo "Supply BAM file as option 2"
	echo "Supply fasta reference genome used to make BAM as option 3"
	echo "Supply bed of TE family locations as option 4"
	echo "Supply TE hierarchy file as option 5"
fi
