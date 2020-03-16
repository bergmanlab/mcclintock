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

	if [ ! -f $outputfolder/../reference/$reference"_locationlist" ]
	then
		mkdir $outputfolder/../reference/$reference
		# Create the individual input files required
		perl splitforRetroSeq.pl $consensus_te_sequences $reference $outputfolder/../reference/
		awk -F".fasta" '{print $1".bed"}' $outputfolder/../reference/$reference"_elementlist" > $outputfolder/../reference/$reference"_locationlist"

		while read element file
		do
			awk -v element=$element '{if ($2==element) print $1}' $te_hierarchy > $outputfolder/../reference/intfile
			while read insert
			do
				grep -w $insert $bed_te_location >> $file
			done < $outputfolder/../reference/intfile
		done < $outputfolder/../reference/$reference"_locationlist"
		rm $outputfolder/../reference/intfile
	fi

	sample=`basename $bam_file .bam`
	mkdir $outputfolder/RetroSeq

	# To run RetroSeq with the more computationally expensive Exonerate step uncomment the next line
	# and comment the alternate discovery step below it
	# bin/retroseq.pl -discover -bam $bam_file -eref $outputfolder/$reference"_elementlist" -output $outputfolder/$sample/$sample.discovery
	bin/retroseq.pl -discover -bam $bam_file -refTEs $outputfolder/../reference/$reference"_locationlist" -output $outputfolder/RetroSeq/$sample.discovery

	bin/retroseq.pl -call -bam $bam_file -input $outputfolder/RetroSeq/$sample.discovery -filter $outputfolder/../reference/$reference"_locationlist" -ref $reference_genome -output $outputfolder/RetroSeq/$sample.calling -orientate yes

	# Extract the relevant results
	echo -e "track name=\"$sample"_RetroSeq"\" description=\"$sample"_RetroSeq"\"" > $outputfolder/RetroSeq/$sample"_retroseq_raw.bed"
	echo -e "track name=\"$sample"_RetroSeq"\" description=\"$sample"_RetroSeq"\"" > $outputfolder/RetroSeq/$sample"_retroseq_redundant.bed"
	echo -e "track name=\"$sample"_RetroSeq"\" description=\"$sample"_RetroSeq"\"" > $outputfolder/RetroSeq/$sample"_retroseq_nonredundant.bed"

	awk '$1!~/#/{print $0}' $outputfolder/RetroSeq/$sample.calling.PE.vcf > $outputfolder/RetroSeq/tmp
	awk -F'[=,\t:]' -v samplename=$sample '{print $1"\t"$11-1"\t"$12"\t"$7"\t"$10"_non-reference_"samplename"_retroseq_rp_\t0\t.\t"$21}' $outputfolder/RetroSeq/tmp > $outputfolder/RetroSeq/$sample"_retroseq_presort.txt"
	bedtools sort -i $outputfolder/RetroSeq/$sample"_retroseq_presort.txt" > $outputfolder/RetroSeq/$sample"_retroseq_sorted_raw.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8"\t"$9}' $outputfolder/RetroSeq/$sample"_retroseq_sorted_raw.txt" > $outputfolder/RetroSeq/$sample"_retroseq_sorted_counted_raw.txt"
	cut -f1-3,5-7 $outputfolder/RetroSeq/$sample"_retroseq_sorted_counted_raw.txt" >> $outputfolder/RetroSeq/$sample"_retroseq_raw.bed"

	# Filter results
	awk '{if ($8 >= 6) print $0}' $outputfolder/RetroSeq/$sample"_retroseq_sorted_counted_raw.txt" > $outputfolder/RetroSeq/$sample"_retroseq_redundant.txt"
	cut -f1-3,5-7 $outputfolder/RetroSeq/$sample"_retroseq_redundant.txt" >> $outputfolder/RetroSeq/$sample"_retroseq_redundant.bed"

	# Filter for redundant predictions
	sort -k1,3 -k4rn $outputfolder/RetroSeq/$sample"_retroseq_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $outputfolder/RetroSeq/tmp
	bedtools sort -i $outputfolder/RetroSeq/tmp >> $outputfolder/RetroSeq/$sample"_retroseq_nonredundant.bed"

	rm $outputfolder/RetroSeq/tmp $outputfolder/RetroSeq/$sample"_retroseq_redundant.txt" $outputfolder/RetroSeq/$sample"_retroseq_sorted_counted_raw.txt" $outputfolder/RetroSeq/$sample"_retroseq_sorted_raw.txt" $outputfolder/RetroSeq/$sample"_retroseq_presort.txt"

else
	echo "Supply TE database as option 1"
	echo "Supply BAM file as option 2"
	echo "Supply fasta reference genome used to make BAM as option 3"
	echo "Supply bed of TE family locations as option 4"
	echo "Supply TE hierarchy file as option 5"
	echo "Supply output folder as option 6"
fi
