#!/usr/bin/bash 

if (( $# > 0 ))
then
	reference=${2##*/}
	reference=${reference%%.*}
	sample=$4
	# Run the relocaTE pipeline
	if [ ! -f $reference"annotation" ]; then
		awk -F'[\t]' '{print $3"\t"$1":"$4".."$5}' $5 > $reference"annotation"
	fi
	perl scripts/relocaTE.pl -t $1 -d $3 -g $2 -1 _1 -2 _2 -o $sample -r $reference"annotation"
	

	# Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$sample"_RelocaTE"\" description=\"$sample"_RelocaTE"\"" > $sample/$sample"_relocate_redundant.bed"
	echo -e "track name=\"$sample"_RelocaTE"\" description=\"$sample"_RelocaTE"\"" > $sample/$sample"_relocate_nonredundant.bed"

	# Extract the relevant information from the output files for each TE and collate them.
	for file in $sample/*/results/*.gff
	do
		awk -F'[\t=;]' -v sample=$sample '$14~/Shared/{print $1"\t"$4"\t"$5"\t"$16+$18"\t"$12"_old_"sample"_relocate_sr_\t0\t."}' $file >> $sample/$sample"_relocate_presort_redundant.txt"
		awk -F'[\t=;.]' -v sample=$sample '$18~/Non-reference/{print $1"\t"$5-1"\t"$6"\t"$20+$22"\t"$13"_new_"sample"_relocate_sr_\t0\t"$9}' $file >> $sample/$sample"_relocate_presort_redundant.txt"
	done

	bedtools sort -i $sample/$sample"_relocate_presort_redundant.txt" > $sample/$sample"_relocate_precut_redundant.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' $sample/$sample"_relocate_precut_redundant.txt" > $sample/$sample"_relocate_counted_precut_redundant.txt"

	cut -f1-3,5-7 $sample/$sample"_relocate_counted_precut_redundant.txt" >> $sample/$sample"_relocate_redundant.bed"

	# Filtering out redundant insertions
    sort -k1,3 -k4rn $sample/$sample"_relocate_counted_precut_redundant.txt" | sort -u -k1,3 | cut -f1-3,5- > $sample/tmp
    bedtools sort -i $sample/tmp >> $sample/$sample"_relocate_nonredundant.bed"

	# Clean up intermediate files
	rm $sample/$sample"_relocate_counted_precut_redundant.txt" $sample/$sample"_relocate_precut_redundant.txt" $sample/$sample"_relocate_presort_redundant.txt" $sample/tmp
	
else
	echo "Supply TE sequence with TSD information in description (format 'TSD=....') as option 1"
	echo "Supply fasta reference file as option 2"
	echo "Supply a directory containing fastq files as option 3"
	echo "Supply an output directory name as option 4"
	echo "Supply reference insertion locations in gff format as option 5"
fi
