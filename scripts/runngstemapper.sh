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

	# Set up base directory for project
	mkdir $outputfolder
	mkdir $outputfolder/$sample/
	mkdir $outputfolder/$sample/samples
	mkdir $outputfolder/$sample/samples/fasta
	mkdir $outputfolder/$sample/reference
	mkdir $outputfolder/$sample/reference/te
	mkdir $outputfolder/$sample/reference/genome
	mkdir $outputfolder/$sample/analysis/
	mkdir $outputfolder/$sample/analysis/align_te
	mkdir $outputfolder/$sample/analysis/fasta_aligned_te
	mkdir $outputfolder/$sample/analysis/align_genome
	mkdir $outputfolder/$sample/analysis/bed_tsd
	mkdir $outputfolder/$sample/analysis/metadata
	mkdir $outputfolder/$sample/analysis/r_data_files
	
	base2=`basename $fasta2_file`
	fasta2=${base2%.*}
	
	cp $consensus_te_seqs $outputfolder/$sample/reference/te
	cp $reference_genome $outputfolder/$sample/reference/genome

	# Run ngs_te_mapper
	R --no-save < sourceCode/ngs_te_mapper.R "$fasta1_file;$fasta2_file" $outputfolder/$sample 1 20

	# Extract only the relevant data from the output file and sort the results
	# Name and description for use with the UCSC genome browser are added to output here.
	awk -F'[\t;]' -v sample=$sample '{print $1"\t"$2"\t"$3"\t"$6"_"$9"_"sample"_ngs_te_mapper_sr_"NR"\t0\t"$5}' $outputfolder/$sample/analysis/bed_tsd/$sample"_1_"$fasta2"insertions.bed" > $outputfolder/$sample/$sample"_ngs_te_mapper_presort.bed"
	echo -e "track name=\"$sample"_ngs_te_mapper"\" description=\"$sample"_ngs_te_mapper"\"" > $outputfolder/$sample/$sample"_ngs_te_mapper_nonredundant.bed"
	bedtools sort -i $outputfolder/$sample/$sample"_ngs_te_mapper_presort.bed" >> $outputfolder/$sample/$sample"_ngs_te_mapper.bed.tmp"
    sed 's/NA/./g' $outputfolder/$sample/$sample"_ngs_te_mapper.bed.tmp" >> $outputfolder/$sample/$sample"_ngs_te_mapper_nonredundant.bed"
	rm $outputfolder/$sample/$sample"_ngs_te_mapper_presort.bed" $outputfolder/$sample/$sample"_ngs_te_mapper.bed.tmp"
	
else
	echo "Supply TE database as option 1"
	echo "Supply Reference genome as option 2"
	echo "Supply sample name as option 3"
	echo "Supply sequence fastas as options 4 and 5"
	echo "Supply output folder as option 6"
	echo "Result will be saved under SAMPLENAME_ngs_te_mapper.bed in the output directory"
fi

