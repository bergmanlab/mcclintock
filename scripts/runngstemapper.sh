#!/bin/bash

if (( $# > 0 ))
then

    # Establish variables
    consensus_te_seqs=$1
    reference_genome=$2
    sample_name=$3
    fasta1_file=$4
    fasta2_file=$5

	# Set up base directory for project
	mkdir $sample_name
	mkdir $sample_name/samples
	mkdir $sample_name/samples/fasta
	mkdir $sample_name/reference
	mkdir $sample_name/reference/te
	mkdir $sample_name/reference/genome
	mkdir $sample_name/analysis/
	mkdir $sample_name/analysis/align_te
	mkdir $sample_name/analysis/fasta_aligned_te
	mkdir $sample_name/analysis/align_genome
	mkdir $sample_name/analysis/bed_tsd
	mkdir $sample_name/analysis/metadata
	mkdir $sample_name/analysis/r_data_files
	
	base2=`basename $fasta2_file`
	fasta2=${base2%.*}
	
	cp $consensus_te_seqs $sample_name/reference/te
	cp $reference_genome $sample_name/reference/genome

	# Run ngs_te_mapper
	R --no-save < sourceCode/ngs_te_mapper.R "$fasta1_file;$fasta2_file" $sample_name 1 20

	# Extract only the relevant data from the output file and sort the results
	# Name and description for use with the UCSC genome browser are added to output here.
	awk -F'[\t;]' -v sample=$3 '{print $1"\t"$2"\t"$3"\t"$6"_"$9"_"sample"_ngs_te_mapper_sr_"NR"\t0\t"$5}' $sample_name/analysis/bed_tsd/$sample_name"_"$fasta2"insertions.bed" > $sample_name/$3"_ngs_te_mapper_presort.bed"
	echo -e "track name=\"$3"_ngs_te_mapper"\" description=\"$3"_ngs_te_mapper"\"" > $sample_name/$3"_ngs_te_mapper_nonredundant.bed"
	bedtools sort -i $sample_name/$3"_ngs_te_mapper_presort.bed" >> $sample_name/$3"_ngs_te_mapper.bed.tmp"
    sed 's/NA/./g' $sample_name/$3"_ngs_te_mapper.bed.tmp" >> $sample_name/$3"_ngs_te_mapper_nonredundant.bed"
	rm $sample_name/$3"_ngs_te_mapper_presort.bed" $sample_name/$3"_ngs_te_mapper.bed.tmp"
	
else
	echo "Supply TE database as option 1"
	echo "Supply Reference genome as option 2"
	echo "Supply sample name as option 3"
	echo "Supply sequence fastas as options 4 and 5"
	echo "Result will be saved under SAMPLENAME_ngs_te_mapper.bed in the output directory"

fi

