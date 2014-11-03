#!/bin/bash

if (( $# > 0 ))
then
	# Set up base directory for project 
	projectdir=$3
	mkdir $projectdir
	mkdir $projectdir/samples
	mkdir $projectdir/samples/fasta
	mkdir $projectdir/reference
	mkdir $projectdir/reference/te
	mkdir $projectdir/reference/genome
	mkdir $projectdir/analysis/
	mkdir $projectdir/analysis/align_te
	mkdir $projectdir/analysis/fasta_aligned_te
	mkdir $projectdir/analysis/align_genome
	mkdir $projectdir/analysis/bed_tsd
	mkdir $projectdir/analysis/metadata
	mkdir $projectdir/analysis/r_data_files
	
	base1=`basename $4`
	samplename=${base1%.*}
	base2=`basename $5`
	fasta2=${base2%.*}
	
	cp $1 $projectdir/reference/te
	cp $2 $projectdir/reference/genome

	# Run ngs_te_mapper
	R --no-save < sourceCode/ngs_te_mapper.R "$4;$5" $projectdir 1 20

	# Extract only the relevant data from the output file and sort the results
	# Name and description for use with the UCSC genome browser are added to output here.
	awk -F'[\t;]' -v sample=$3 '{print $1"\t"$2"\t"$3"\t"$6"_"$9"_"sample"_ngs_te_mapper_sr\t0\t"$5}' $projectdir/analysis/bed_tsd/$samplename"_"$fasta2"insertions.bed" > $projectdir/$3"_ngs_te_mapper_presort.bed"
	echo -e "track name=\"$3"_ngs_te_mapper"\" description=\"$3"_ngs_te_mapper"\"" > $projectdir/$3"_ngs_te_mapper.bed"
	bedtools sort -i $projectdir/$3"_ngs_te_mapper_presort.bed" >> $projectdir/$3"_ngs_te_mapper.bed.tmp"
    sed 's/NA/./g' $projectdir/$3"_ngs_te_mapper.bed.tmp" >> $projectdir/$3"_ngs_te_mapper.bed"
	rm $projectdir/$3"_ngs_te_mapper_presort.bed" $projectdir/$3"_ngs_te_mapper.bed.tmp"
	
else
	echo "Supply TE database as option 1"
	echo "Supply Reference genome as option 2"
	echo "Supply output directory as option 3"
	echo "Supply sequence fastas as options 4 and 5"
	echo "Result will be saved under SAMPLENAME_ngs_te_mapper.bed in the output directory"

fi

