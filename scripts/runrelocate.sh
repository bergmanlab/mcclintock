#!/usr/bin/bash 

if (( $# > 0 ))
then
	# Run the relocaTE pipeline
	awk -F'[\t=]' '{print $9"\t"$1":"$4".."$5}' $5 > annotation
	perl scripts/relocaTE.pl -t $1 -d $3 -g $2 -1 _1 -2 _2 -o $4 -r annotation
	# Extract the relevant information from the output files for each TE and collate them.
	# Name and description for use with the UCSC genome browser are added to output here.
        echo -e "track name=\"RelocaTE\" description=\"RelocaTE\"" >> $4/relocate_results.bed
	for file in $4/*/results/*.gff 
	do
		awk -F'[\t=;]' '$14~/Shared/{print $1"\t"$4"\t"$5"\t"$12"_old\t0\t."}' $file >> $4/relocate_results.bed
		awk -F'[\t=;.]' '$18~/Non-reference/{print $1"\t"$5"\t"$6"\t"$13"_new\t0\t"$9}' $file >> $4/relocate_results.bed
	done	
	
else
	echo "Supply TE sequence with TSD information in description (format 'TSD=....') as option 1"
	echo "Supply fasta reference file as option 2"
	echo "Supply a directory containing fastq files as option 3"
	echo "Supply an output directory name as option 4"
	echo "Supply reference insertion locations in gff format as option 5"
fi
