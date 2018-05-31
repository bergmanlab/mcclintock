#!/bin/bash

if (( $# > 0 ))
then

	# Establish variables
	consensus_te_seqs=$1
	reference_genome=$2
	sample=$3
	fasta1_file=$4
	fasta2_file=$5
	outputfolder=$6
	test_dir=`pwd`
	
	mkdir -p $outputfolder/ngs_te_mapper

	if [[ $fasta2_file == "false" ]]
	then
		base1=`basename $fasta1_file`
		samplename=${base1%%.*}
		# Run ngs_te_mapper on single end data
		Rscript --vanilla $test_dir/sourceCode/ngs_te_mapper.R sample=$fasta1_file genome=$reference_genome teFile=$consensus_te_seqs tsd=20 output=$outputfolder/ngs_te_mapper sourceCodeFolder=$test_dir/sourceCode
	else
		# Run ngs_te_mapper with both fastq files
		Rscript --vanilla $test_dir/sourceCode/ngs_te_mapper.R sample=$fasta1_file\;$fasta2_file genome=$reference_genome teFile=$consensus_te_seqs tsd=20 output=$outputfolder/ngs_te_mapper sourceCodeFolder=$test_dir/sourceCode
	fi

	# Extract only the relevant data from the output file and sort the results
	# Name and description for use with the UCSC genome browser are added to output here.
	if [[ $fasta2_file == "false" ]]
	then
		awk -F'[\t;]' -v sample=$sample '{print $1"\t"$2"\t"$3"\t"$6"_"$9"_"sample"_ngs_te_mapper_sr_"NR"\t0\t"$5}' $outputfolder/ngs_te_mapper/bed_tsd/$samplename"insertions.bed" > $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper_presort.bed"
	else
		base1=`basename $fasta1_file`
		samplename1=${base1%%_1.f*}
		samplename1=${samplename1%%.*}
		base2=`basename $fasta2_file`
		samplename2=${base2%%_1.f*}
		samplename2=${samplename2%%.*}
		awk -F'[\t;]' -v sample=$sample '{print $1"\t"$2"\t"$3"\t"$6"_"$9"_"sample"_ngs_te_mapper_sr_"NR"\t0\t"$5}' $outputfolder/ngs_te_mapper/bed_tsd/$samplename1"_"$samplename2"insertions.bed" > $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper_presort.bed"
	fi

	echo -e "track name=\"$sample"_ngs_te_mapper"\" description=\"$sample"_ngs_te_mapper"\"" > $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper_nonredundant.bed"
	bedtools sort -i $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper_presort.bed" >> $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper.bed.tmp"
	sed 's/NA/./g' $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper.bed.tmp" >> $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper_nonredundant.bed"
	rm $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper_presort.bed" $outputfolder/ngs_te_mapper/$sample"_ngs_te_mapper.bed.tmp"
	
else
	echo "Supply TE database as option 1"
	echo "Supply Reference genome as option 2"
	echo "Supply sample name as option 3"
	echo "Supply sequence fastas as options 4 and 5"
	echo "Supply output folder as option 6"
	echo "Result will be saved under SAMPLENAME_ngs_te_mapper.bed in the output directory"
fi

