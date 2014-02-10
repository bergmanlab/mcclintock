#!/bin/bash

if (( $# > 0 ))
then
	#set up base directory for project 
	projectdir=$3
	mkdir $projectdir
	mkdir $projectdir/samples
	mkdir $projectdir/samples/fasta
	mkdir $projectdir/reference
	mkdir $projectdir/reference/te
	mkdir $projectdir/reference/genome
	mkdir $projectdir/analysis/
	mkdir $projectdir/analysis/psl_te
	mkdir $projectdir/analysis/fasta_aligned_te
	mkdir $projectdir/analysis/psl_genome
	mkdir $projectdir/analysis/bed_tsd
	mkdir $projectdir/analysis/metadata
	mkdir $projectdir/analysis/r_data_files
	mkdir $projectdir/logo
	
	base1=`basename $4`
	samplename=${base1%.*}
	base2=`basename $5`
	fasta2=${base2%.*}
	
	cp $1 $projectdir/reference/te
	cp $2 $projectdir/reference/genome

	#run ngs_te_mapper on different files has if it was only one sample (for paired end)
	#the names of the files have to be separated by ";"
	R --no-save < sourceCode/ngs_te_mapper.R "$4;$5" $projectdir 1 20

	R --no-save < sourceCode/ngs_te_logo.R $projectdir 25
	#Extract only the relevant data from the output file
	#Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"ngs_te_mapper\" description=\"ngs_te_mapper\"" >> $projectdir/$samplename"_ngs_te_mapper.bed"
	awk -F'[\t;]' '{print $1"\t"$2"\t"$3"\t"$9"_new\t0\t"$8}' $projectdir/analysis/bed_tsd/$samplename"_"$fasta2"insertions.bed" >> $projectdir/$samplename"_ngs_te_mapper.bed"
else
	echo "Supply TE database as option 1"
	echo "Supply Reference genome as option 2"
	echo "Supply output directory as option 3"
	echo "Supply sequence fastas as options 4 and 5"
	echo "Result will be saved under SAMPLENAME_ngs_te_mapper.bed in the output directory"
fi

