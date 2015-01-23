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
	outputfolder=$6

	if [ ! -d $outputfolder/$reference ]
	then
		mkdir -p $outputfolder
		mkdir $outputfolder/$reference
		# Create the individual input files required
		perl splitforRetroSeq.pl $consensus_te_sequences $reference $outputfolder
		awk -F".fa" '{print $1".bed"}' $outputfolder/$reference"_elementlist" > $outputfolder/$reference"_locationlist"

		while read element file
		do
			awk -v element=$element '{if ($2==element) print $1}' $te_hierarchy > $outputfolder/intfile
			while read insert
			do
				grep -w $insert $bed_te_location >> $file
			done < $outputfolder/intfile
		done < $outputfolder/$reference"_locationlist"
		rm $outputfolder/intfile
	fi

	sample=`basename $bam_file .bam`
	mkdir $outputfolder/$sample

	# To run RetroSeq with the more computationally expensive Exonerate step uncomment the next line
	# and comment the alternate discovery step below it
	# bin/retroseq.pl -discover -bam $bam_file -eref $outputfolder/$reference"_elementlist" -output $outputfolder/$sample/$sample.discovery
	bin/retroseq.pl -discover -bam $bam_file -refTEs $outputfolder/$reference"_locationlist" -output $outputfolder/$sample/$sample.discovery

	bin/retroseq.pl -call -bam $bam_file -input $outputfolder/$sample/$sample.discovery -filter $outputfolder/$reference"_locationlist" -ref $reference_genome -output $outputfolder/$sample/$sample.calling -orientate yes

	# Extract the relevant results
	echo -e "track name=\"$sample"_RetroSeq"\" description=\"$sample"_RetroSeq"\"" > $outputfolder/$sample/$sample"_retroseq_raw.bed"
	echo -e "track name=\"$sample"_RetroSeq"\" description=\"$sample"_RetroSeq"\"" > $outputfolder/$sample/$sample"_retroseq_redundant.bed"
	echo -e "track name=\"$sample"_RetroSeq"\" description=\"$sample"_RetroSeq"\"" > $outputfolder/$sample/$sample"_retroseq_nonredundant.bed"

	awk '$1!~/#/{print $0}' $outputfolder/$sample/$sample.calling.PE.vcf > $outputfolder/$sample/tmp
	awk -F'[=,\t:]' -v samplename=$sample '{print $1"\t"$11"\t"$12"\t"$6"\t"$10"_new_"samplename"_retroseq_rp_\t0\t.\t"$21}' $outputfolder/$sample/tmp > $outputfolder/$sample/$sample"_retroseq_presort.txt"
	bedtools sort -i $outputfolder/$sample/$sample"_retroseq_presort.txt" > $outputfolder/$sample/$sample"_retroseq_sorted_raw.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8}' $outputfolder/$sample/$sample"_retroseq_sorted_raw.txt" > $outputfolder/$sample/$sample"_retroseq_sorted_counted_raw.txt"
	cut -f1-3,5-7 $outputfolder/$sample/$sample"_retroseq_sorted_counted_raw.txt" >> $outputfolder/$sample/$sample"_retroseq_raw.bed"

	# Filter results
	awk '{if ($8 >= 6) print $0}' $outputfolder/$sample/$sample"_retroseq_sorted_counted_raw.txt" > $outputfolder/$sample/$sample"_retroseq_redundant.txt"
	cut -f1-3,5-7 $outputfolder/$sample/$sample"_retroseq_redundant.txt" >> $outputfolder/$sample/$sample"_retroseq_redundant.bed"

	# Filter for redundant predictions
	sort -k1,3 -k4rn $outputfolder/$sample/$sample"_retroseq_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $outputfolder/$sample/tmp
	bedtools sort -i $outputfolder/$sample/tmp >> $outputfolder/$sample/$sample"_retroseq_nonredundant.bed"

	rm $outputfolder/$sample/tmp $outputfolder/$sample/$sample"_retroseq_redundant.txt" $outputfolder/$sample/$sample"_retroseq_sorted_counted_raw.txt" $outputfolder/$sample/$sample"_retroseq_sorted_raw.txt" $outputfolder/$sample/$sample"_retroseq_presort.txt"

else
	echo "Supply TE database as option 1"
	echo "Supply BAM file as option 2"
	echo "Supply fasta reference genome used to make BAM as option 3"
	echo "Supply bed of TE family locations as option 4"
	echo "Supply TE hierarchy file as option 5"
	echo "Supply output folder as option 6"
fi
