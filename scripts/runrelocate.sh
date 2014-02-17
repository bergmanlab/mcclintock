#!/usr/bin/bash 

if (( $# > 0 ))
then
	# Run the relocaTE pipeline
	awk -F'[\t]' '{print $3"\t"$1":"$4".."$5}' $5 > annotation
	perl scripts/relocaTE.pl -t $1 -d $3 -g $2 -1 _1 -2 _2 -o $4 -r annotation
	
	# Extract the relevant information from the output files for each TE and collate them.
	# Name and description for use with the UCSC genome browser are added to output here.
	for file in $4/*/results/*.gff 
	do
		awk -F'[\t=;]' '$14~/Shared/{print $1"\t"$4"\t"$5"\t"$12"_old\t0\t."}' $file >> $4/$4_relocate_presort.bed
		awk -F'[\t=;.]' '$18~/Non-reference/{print $1"\t"$5"\t"$6"\t"$13"_new\t0\t"$9}' $file >> $4/$4_relocate_presort.bed
	done	
	echo -e "track name=\"$4"_RelocaTE"\" description=\"$4"_RelocaTE"\"" > $4/$4_relocate.bed
	bedtools sort -i $4/$4_relocate_presort.bed >> $4/$4_relocate.bed
	rm $4/$4_relocate_presort.bed
	
else
	echo "Supply TE sequence with TSD information in description (format 'TSD=....') as option 1"
	echo "Supply fasta reference file as option 2"
	echo "Supply a directory containing fastq files as option 3"
	echo "Supply an output directory name as option 4"
	echo "Supply reference insertion locations in gff format as option 5"
fi
