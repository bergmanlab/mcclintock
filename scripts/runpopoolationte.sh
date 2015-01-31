#!/usr/bin/bash

if (( $# > 0 ))
then

	# Establish variables
	reference_genome=$1
	te_hierarchy=$2
	fastq1=$3
	fastq2=$4
	gff_te_locations=$5
	processors=$6
	outputfolder=$7

	# Get the name of the sample being analysed (catch either fq or fastq extensions).
	samplename=`basename $fastq1 _1.fastq`
	sample=`basename $samplename _1.fq`

	# Create a directory for analysis files
	mkdir -p $outputfolder
	mkdir $outputfolder/$sample

	# Change the labels of reads in the fastq files to ensure there are no spaces and they end in /1 or /2
	awk -v samplen=$sample '{if(NR%4==1) $0=sprintf("@"samplen".%d/1",(1+i++)); print;}' $fastq1 > $outputfolder/$sample/reads1.fastq
	awk -v samplen=$sample '{if(NR%4==1) $0=sprintf("@"samplen".%d/2",(1+i++)); print;}' $fastq2 > $outputfolder/$sample/reads2.fastq

	# Perform 2 bwa alignments.
	bwa bwasw -t $processors $reference_genome $outputfolder/$sample/reads1.fastq > $outputfolder/$sample/1.sam
	bwa bwasw -t $processors $reference_genome $outputfolder/$sample/reads2.fastq > $outputfolder/$sample/2.sam

	# Combine each alignment using included script.
	perl samro.pl --sam1 $outputfolder/$sample/1.sam --sam2 $outputfolder/$sample/2.sam --fq1 $outputfolder/$sample/reads1.fastq --fq2 $outputfolder/$sample/reads2.fastq --output $outputfolder/$sample/pe-reads.sam

	rm $outputfolder/$sample/reads1.fastq
	rm $outputfolder/$sample/reads2.fastq

	# Sort the alignment and remove intermediate files to save space.
	samtools view -t $reference_genome".fai" -Sb $outputfolder/$sample/pe-reads.sam | samtools sort - $outputfolder/$sample/pe-reads.sorted
	rm $outputfolder/$sample/pe-reads.sam
	samtools view $outputfolder/$sample/pe-reads.sorted.bam > $outputfolder/$sample/pe-reads.sorted.sam

	# Perform the rest of the PoPoolationTE workflow
	perl identify-te-insertsites.pl --input $outputfolder/$sample/pe-reads.sorted.sam --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --narrow-range 75 --min-count 3 --min-map-qual 15 --output $outputfolder/$sample/te-fwd-rev.txt
	perl genomic-N-2gtf.pl --input $reference_genome > $outputfolder/$sample/poly_n.gtf
	perl crosslink-te-sites.pl --directional-insertions $outputfolder/$sample/te-fwd-rev.txt --min-dist 74 --max-dist 250 --output $outputfolder/$sample/te-inserts.txt --single-site-shift 100 --poly-n $outputfolder/$sample/poly_n.gtf --te-hierarchy $te_hierarchy --te-hier-level family

	# Convert the gff format to a format that PoPoolationTE understads.
	awk -F'[\t=;]' '{if($7=="-")print $1"\tF\t"$5"\t"$10"\n"$1"\tR\t"$4"\t"$10; else print $1"\tF\t"$4"\t"$10"\n"$1"\tR\t"$5"\t"$10}' $gff_te_locations > $outputfolder/$sample/known-te-insertions.txt
	perl update-teinserts-with-knowntes.pl --known $outputfolder/$sample/known-te-insertions.txt --output $outputfolder/$sample/te-insertions.updated.txt --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --max-dist 300 --te-insertions $outputfolder/$sample/te-inserts.txt --single-site-shift 100
	perl estimate-polymorphism.pl --sam-file $outputfolder/$sample/pe-reads.sorted.sam --te-insert-file $outputfolder/$sample/te-insertions.updated.txt --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --min-map-qual 15 --output $outputfolder/$sample/te-polymorphism.txt
	perl filter-teinserts.pl --te-insertions $outputfolder/$sample/te-polymorphism.txt --output $outputfolder/$sample/te-poly-filtered.txt --discard-overlapping --min-count 10

	# Create an output file in BED format selecting only the relevant inserts.
	# Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$sample"_PoPoolationTE"\" description=\"$sample"_PoPoolationTE"\"" > $outputfolder/$sample/$sample"_popoolationte_raw.bed"
	echo -e "track name=\"$sample"_PoPoolationTE"\" description=\"$sample"_PoPoolationTE"\"" > $outputfolder/$sample/$sample"_popoolationte_redundant.bed"
	echo -e "track name=\"$sample"_PoPoolationTE"\" description=\"$sample"_PoPoolationTE"\"" > $outputfolder/$sample/$sample"_popoolationte_nonredundant.bed"

	# Convert results to bed with no filtering
	awk -F"\t" -v sample=$sample '{if($7!~/-/) printf "%s\t%d\t%d\t%d\t%s%s%s%s\t0\t.\t%s\n",$1,$2,$2,$13+$20,$4,"_old_",sample,"_popoolationte_rp_",$3; else printf "%s\t%d\t%d\t%d\t%s%s%s%s\t0\t.\t%s\n",$1,$2,$2,$13+$20,$4,"_new_",sample,"_popoolationte_rp_",$3}' $outputfolder/$sample/te-poly-filtered.txt > $outputfolder/$sample/$sample"_popoolationte_presort_raw.txt"
    sort -k1,1 -k2,2n $outputfolder/$sample/$sample"_popoolationte_presort_raw.txt" $outputfolder/$sample/$sample"_popoolationte_precut_raw.txt" > $outputfolder/$sample/$sample"_popoolationte_precut_raw.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8}' $outputfolder/$sample/$sample"_popoolationte_precut_raw.txt" > $outputfolder/$sample/$sample"_popoolationte_counted_precut_raw.txt"
	cut -f1-3,5-7 $outputfolder/$sample/$sample"_popoolationte_counted_precut_raw.txt" >> $outputfolder/$sample/$sample"_popoolationte_raw.bed"

	# Filtering for inserts with support for both ends
	awk '{if($8=="FR") print $0}' $outputfolder/$sample/$sample"_popoolationte_counted_precut_raw.txt" > $outputfolder/$sample/$sample"_popoolationte_counted_precut_redundant.txt"
	cut -f1-3,5-7 $outputfolder/$sample/$sample"_popoolationte_counted_precut_redundant.txt" >> $outputfolder/$sample/$sample"_popoolationte_redundant.bed"

	# Filtering out redundant insertions
	sort -k1,3 -k4rn $outputfolder/$sample/$sample"_popoolationte_counted_precut_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $outputfolder/$sample/tmp
	sort -k1,1 -k2,2n $outputfolder/$sample/tmp >> $outputfolder/$sample/$sample"_popoolationte_nonredundant.bed"

	# Clean up intermediate files
	rm $outputfolder/$sample/tmp $outputfolder/$sample/$sample"_popoolationte_counted_precut_raw.txt" $outputfolder/$sample/$sample"_popoolationte_counted_precut_redundant.txt" $outputfolder/$sample/$sample"_popoolationte_precut_raw.txt" $outputfolder/$sample/$sample"_popoolationte_presort_raw.txt" $outputfolder/$sample/1.sam $outputfolder/$sample/2.sam $outputfolder/$sample/pe-reads.sorted.bam

else
	echo "Supply fasta reference file as option 1"
	echo "Supply TE hierarchy file as option 2"
	echo "Supply fastq1 as option 3"
	echo "Supply fastq2 as option 4"
	echo "Supply a gff3 annotation of the TEs in the genome as option 5"
	echo "Supply a number of processors as option 6"
	echo "Supply output folder as option 7"
fi

