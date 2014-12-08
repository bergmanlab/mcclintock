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

    # Get the name of the sample being analysed.
    samplename=`basename $fastq1 _1.fastq`

    # Create a directory for analysis files
    mkdir $samplename

	# Remove any spaces in the names of reads to prevent errors.
	tr -d ' ' < $fastq1 > $samplename/reads1.fq
	tr -d ' ' < $fastq2 > $samplename/reads2.fq	

	# Perform 2 bwa alignments.
	bwa bwasw -t $processors $reference_genome $samplename/reads1.fq > $samplename/1.sam
	bwa bwasw -t $processors $reference_genome  $samplename/reads2.fq > $samplename/2.sam

	# Combine each alignment using included script.
	perl samro.pl --sam1 $samplename/1.sam --sam2 $samplename/2.sam --fq1 $samplename/reads1.fq --fq2 $samplename/reads2.fq --output $samplename/pe-reads.sam

	rm $samplename/reads1.fq
	rm $samplename/reads2.fq

	# Sort the alignment and remove intermediate files to save space.
	samtools view -Sb $samplename/pe-reads.sam | samtools sort - $samplename/pe-reads.sorted
	rm $samplename/pe-reads.sam
	samtools view $samplename/pe-reads.sorted.bam > $samplename/pe-reads.sorted.sam

	# Perform the rest of the PoPoolationTE workflow
	perl identify-te-insertsites.pl --input $samplename/pe-reads.sorted.sam --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --narrow-range 75 --min-count 3 --min-map-qual 15 --output $samplename/te-fwd-rev.txt
	perl genomic-N-2gtf.pl --input $reference_genome > $samplename/poly_n.gtf
	perl crosslink-te-sites.pl --directional-insertions $samplename/te-fwd-rev.txt --min-dist 74 --max-dist 250 --output $samplename/te-inserts.txt --single-site-shift 100 --poly-n $samplename/poly_n.gtf --te-hierarchy $te_hierarchy --te-hier-level family

	# Convert the gff format to a format that PoPoolationTE understads.
	awk -F'[\t=;]' '{if($7=="-")print $1"\tF\t"$5"\t"$10"\n"$1"\tR\t"$4"\t"$10; else print $1"\tF\t"$4"\t"$10"\n"$1"\tR\t"$5"\t"$10}' $gff_te_locations > $samplename/known-te-insertions.txt
	perl update-teinserts-with-knowntes.pl --known $samplename/known-te-insertions.txt --output $samplename/te-insertions.updated.txt --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --max-dist 300 --te-insertions $samplename/te-inserts.txt --single-site-shift 100
	perl estimate-polymorphism.pl --sam-file $samplename/pe-reads.sorted.sam --te-insert-file $samplename/te-insertions.updated.txt --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --min-map-qual 15 --output $samplename/te-polymorphism.txt
	perl filter-teinserts.pl --te-insertions $samplename/te-polymorphism.txt --output $samplename/te-poly-filtered.txt --discard-overlapping --min-count 10

	# Create an output file in BED format selecting only the relevant inserts.
	# Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$samplename"_PoPoolationTE"\" description=\"$samplename"_PoPoolationTE"\"" > $samplename/$samplename"_popoolationte_raw.bed"
	echo -e "track name=\"$samplename"_PoPoolationTE"\" description=\"$samplename"_PoPoolationTE"\"" > $samplename/$samplename"_popoolationte_redundant.bed"
	echo -e "track name=\"$samplename"_PoPoolationTE"\" description=\"$samplename"_PoPoolationTE"\"" > $samplename/$samplename"_popoolationte_nonredundant.bed"

	# Convert results to bed with no filtering
	awk -F"\t" -v sample=$samplename '{if($7!~/-/) printf "%s\t%d\t%d\t%d\t%s%s%s%s\t0\t.\t%s\n",$1,$2,$2,$13+$20,$4,"_old_",sample,"_popoolationte_rp_",$3; else printf "%s\t%d\t%d\t%d\t%s%s%s%s\t0\t.\t%s\n",$1,$2,$2,$13+$20,$4,"_new_",sample,"_popoolationte_rp_",$3}' $samplename/te-poly-filtered.txt > $samplename/$samplename"_popoolationte_presort_raw.txt"
    sort -k1,1 -k2,2n $samplename/$samplename"_popoolationte_presort_raw.txt" $samplename/$samplename"_popoolationte_precut_raw.txt" > $samplename/$samplename"_popoolationte_precut_raw.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8}' $samplename/$samplename"_popoolationte_precut_raw.txt" > $samplename/$samplename"_popoolationte_counted_precut_raw.txt"
	cut -f1-3,5-7 $samplename/$samplename"_popoolationte_counted_precut_raw.txt" >> $samplename/$samplename"_popoolationte_raw.bed"

	# Filtering for inserts with support for both ends
	awk '{if($8=="FR") print $0}' $samplename/$samplename"_popoolationte_counted_precut_raw.txt" > $samplename/$samplename"_popoolationte_counted_precut_redundant.txt"
	cut -f1-3,5-7 $samplename/$samplename"_popoolationte_counted_precut_redundant.txt" >> $samplename/$samplename"_popoolationte_redundant.bed"

	# Filtering out redundant insertions
	sort -k1,3 -k4rn $samplename/$samplename"_popoolationte_counted_precut_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $samplename/tmp
	sort -k1,1 -k2,2n $samplename/tmp >> $samplename/$samplename"_popoolationte_nonredundant.bed"

	# Clean up intermediate files
	rm $samplename/tmp $samplename/$samplename"_popoolationte_counted_precut_raw.txt" $samplename/$samplename"_popoolationte_counted_precut_redundant.txt" $samplename/$samplename"_popoolationte_precut_raw.txt" $samplename/$samplename"_popoolationte_presort_raw.txt" $samplename/1.sam $samplename/2.sam $samplename/pe-reads.sorted.bam

else
	echo "Supply fasta reference file as option 1"
	echo "Supply TE hierarchy file as option 2"
	echo "Supply fastq1 as option 3"
	echo "Supply fastq2 as option 4"
	echo "Supply a gff3 annotation of the TEs in the genome as option 5"
	echo "Supply a number of processors as option 6"
fi

