#!/usr/bin/bash

if (( $# > 0 ))
then
	# Get the name of the sample being analysed.
	samplename=`basename $4 _1.fastq`
	
	# Create a directory for analysis files
	mkdir $samplename
	
	# Mask any copies of the TE present in the reference genome.
	if [!-f $1.masked.bwt]; then
		RepeatMasker -no_is -nolow -norna --lib $2 -pa 4 $1
		cat $1.masked $2 > combinedref.fa
		mv combinedref.fa $1.masked
		bwa index $1.masked
	fi

	# Remove any spaces in the names of reads to prevent errors.
	tr -d ' ' < $4 > $samplename/reads1.fq
	tr -d ' ' < $5 > $samplename/reads2.fq	

	# Perform 2 bwa alignments.
	bwa bwasw -t 6 $1.masked $samplename/reads1.fq > $samplename/1.sam
	bwa bwasw -t 6 $1.masked $samplename/reads2.fq > $samplename/2.sam

	# Combine each alignment using included script.
	perl samro.pl --sam1 $samplename/1.sam --sam2 $samplename/2.sam --fq1 $samplename/reads1.fq --fq2 $samplename/reads2.fq --output $samplename/pe-reads.sam

	# Sort the alignment and remove intermediate files to save space.
	#rm 1.sam
	#rm 2.sam
	samtools view -Sb $samplename/pe-reads.sam | samtools sort - $samplename/pe-reads.sorted
	#rm pe-reads.sam
	samtools view $samplename/pe-reads.sorted.bam > $samplename/pe-reads.sorted.sam
	#rm pe-reads.bam

	# Perform the rest of the PoPoolationTE workflow
	perl identify-te-insertsites.pl --input $samplename/pe-reads.sorted.sam --te-hierarchy-file $3 --te-hierarchy-level family --narrow-range 75 --min-count 3 --min-map-qual 15 --output $samplename/te-fwd-rev.txt
	perl genomic-N-2gtf.pl --input $1.masked > $samplename/poly_n.gtf
	perl crosslink-te-sites.pl --directional-insertions $samplename/te-fwd-rev.txt --min-dist 74 --max-dist 250 --output $samplename/te-inserts.txt --single-site-shift 100 --poly-n $samplename/poly_n.gtf --te-hierarchy $3 --te-hier-level family

	# Convert the gff format to a format that PoPoolationTE understads.
	awk -F'[\t=;]' '{if($7=="-")print $1"\tF\t"$5"\t"$10"\n"$1"\tR\t"$4"\t"$10; else print $1"\tF\t"$4"\t"$10"\n"$1"\tR\t"$5"\t"$10}' $6 > $samplename/known-te-insertions.txt
	perl update-teinserts-with-knowntes.pl --known $samplename/known-te-insertions.txt --output $samplename/te-insertions.updated.txt --te-hierarchy-file $3 --te-hierarchy-level family --max-dist 300 --te-insertions $samplename/te-inserts.txt --single-site-shift 100
	perl estimate-polymorphism.pl --sam-file $samplename/pe-reads.sorted.sam --te-insert-file $samplename/te-insertions.updated.txt --te-hierarchy-file $3 --te-hierarchy-level family --min-map-qual 15 --output $samplename/te-polymorphism.txt

	# Create an output file in BED format selecting only the relevant inserts.
	# Name and description for use with the UCSC genome browser are added to output here.
	awk -F"\t" '{if($3=="FR" && $7!~/-/)print $1"\t"$2"\t"$2"\t"$4"_old\t0\t."; else if($3=="FR") print $1"\t"$2"\t"$2"\t"$4"_new\t0\t."}' $samplename/te-polymorphism.txt >> $samplename/$samplename"_popoolationte_presort.bed"
	echo -e "track name=\"$samplename"_PoPoolationTE"\" description=\"$samplename"_PoPoolationTE"\"" > $samplename/$samplename"_popoolationte.bed"
	sort -k1,1 -k2,2n $samplename/$samplename"_popoolationte_presort.bed" >> $samplename/$samplename"_popoolationte.bed"
	rm $samplename/$samplename"_popoolationte_presort.bed"

else
	echo "Supply fasta reference file as option 1"
	echo "Supply fasta TE database file as option 2"
	echo "Supply TE hierarchy file as option 3"
	echo "Supply fastq1 as option 4"
	echo "Supply fastq2 as option 5" 
	echo "Supply a gff3 annotation of the TEs in the genome as option 6"
fi

