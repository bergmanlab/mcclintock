#!/usr/bin/bash 

if (( $# > 0 ))
then

	# Establish variables
	relocate_te_sequences=$1
	reference_genome=$2
	reference=${2##*/}
	reference=${reference%%.*}
	fastq_directory=$3
	sample=$4
	gff_te_locations=$5
	outputfolder=$6
	test_dir=`pwd`

	mkdir -p $outputfolder

	# Run the relocaTE pipeline
	if [ ! -f $outputfolder/$reference"annotation" ]; then
		awk -F'[\t]' '{print $3"\t"$1":"$4".."$5}' $gff_te_locations > $outputfolder/$reference"annotation"
	fi
	cd $outputfolder
	perl $test_dir/scripts/relocaTE.pl -t $relocate_te_sequences -d $fastq_directory -g $reference_genome -1 _1 -2 _2 -o $sample -r $reference"annotation"
	cd ..

	# Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$sample"_RelocaTE"\" description=\"$sample"_RelocaTE"\"" > $outputfolder/$sample/$sample"_relocate_redundant.bed"
	echo -e "track name=\"$sample"_RelocaTE"\" description=\"$sample"_RelocaTE"\"" > $outputfolder/$sample/$sample"_relocate_nonredundant.bed"

	# Extract the relevant information from the output files for each TE and collate them.
	for file in $outputfolder/$sample/*/results/*.gff
	do
		awk -F'[\t=;]' -v sample=$sample '$14~/Shared/{print $1"\t"$4"\t"$5"\t"$16+$18"\t"$12"_old_"sample"_relocate_sr_\t0\t."}' $file >> $outputfolder/$sample/$sample"_relocate_presort_redundant.txt"
		awk -F'[\t=;.]' -v sample=$sample '$18~/Non-reference/{print $1"\t"$5-1"\t"$6"\t"$20+$22"\t"$13"_new_"sample"_relocate_sr_\t0\t"$9}' $file >> $outputfolder/$sample/$sample"_relocate_presort_redundant.txt"
	done

	bedtools sort -i $outputfolder/$sample/$sample"_relocate_presort_redundant.txt" > $outputfolder/$sample/$sample"_relocate_precut_redundant.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7}' $outputfolder/$sample/$sample"_relocate_precut_redundant.txt" > $outputfolder/$sample/$sample"_relocate_counted_precut_redundant.txt"

	cut -f1-3,5-7 $outputfolder/$sample/$sample"_relocate_counted_precut_redundant.txt" >> $outputfolder/$sample/$sample"_relocate_redundant.bed"

	# Filtering out redundant insertions
	sort -k1,3 -k4rn $outputfolder/$sample/$sample"_relocate_counted_precut_redundant.txt" | sort -u -k1,3 | cut -f1-3,5- > $outputfolder/$sample/tmp
	bedtools sort -i $outputfolder/$sample/tmp >> $outputfolder/$sample/$sample"_relocate_nonredundant.bed"

	# Clean up intermediate files
	rm $outputfolder/$sample/$sample"_relocate_counted_precut_redundant.txt" $outputfolder/$sample/$sample"_relocate_precut_redundant.txt" $outputfolder/$sample/$sample"_relocate_presort_redundant.txt" $outputfolder/$sample/tmp
	
else
	echo "Supply TE sequence with TSD information in description (format 'TSD=....') as option 1"
	echo "Supply fasta reference file as option 2"
	echo "Supply a directory containing fastq files as option 3"
	echo "Supply an output directory name as option 4"
	echo "Supply reference insertion locations in gff format as option 5"
	echo "Supply output folder as option 6"
fi
