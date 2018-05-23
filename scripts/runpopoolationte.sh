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
	samplename=`basename $fastq1`
	sample=${samplename%%.f*}

	# Create a directory for analysis files
	mkdir $outputfolder/PoPoolationTE

	# Change the labels of reads in the fastq files to ensure there are no spaces and they end in /1 or /2
	awk -v samplen=$sample '{if(NR%4==1) $0=sprintf("@"samplen".%d/1",(1+i++)); print;}' $fastq1 > $outputfolder/PoPoolationTE/reads1.fastq
	awk -v samplen=$sample '{if(NR%4==1) $0=sprintf("@"samplen".%d/2",(1+i++)); print;}' $fastq2 > $outputfolder/PoPoolationTE/reads2.fastq

	# Get the read length - get average of read lengths in case of differing lengths
	readlength=`sed '2q;d' $outputfolder/PoPoolationTE/reads1.fastq | wc | awk '{print $3-1}'`
	readlength2=`sed '2q;d' $outputfolder/PoPoolationTE/reads2.fastq | wc | awk '{print $3-1}'`
	readlength=$((($readlength + $readlength2) / 2))

	# Perform 2 bwa alignments.
	bwa bwasw -t $processors $reference_genome $outputfolder/PoPoolationTE/reads1.fastq > $outputfolder/PoPoolationTE/1.sam
	bwa bwasw -t $processors $reference_genome $outputfolder/PoPoolationTE/reads2.fastq > $outputfolder/PoPoolationTE/2.sam

	# Combine each alignment using included script.
	perl samro.pl --sam1 $outputfolder/PoPoolationTE/1.sam --sam2 $outputfolder/PoPoolationTE/2.sam --fq1 $outputfolder/PoPoolationTE/reads1.fastq --fq2 $outputfolder/PoPoolationTE/reads2.fastq --output $outputfolder/PoPoolationTE/pe-reads.sam

	median_insertsize=`cut -f9 $outputfolder/PoPoolationTE/pe-reads.sam | sort -S $memory"G" -n | awk '{if ($1 > 0) ins[reads++]=$1; } END { print ins[int(reads/2)]; }'`
	printf "\nMedian insert size = $median_insertsize\n\n" | tee -a /dev/stderr

	maxdist="$(($median_insertsize * 3 + $readlength))"

	rm $outputfolder/PoPoolationTE/reads1.fastq
	rm $outputfolder/PoPoolationTE/reads2.fastq

	# Sort the alignment and remove intermediate files to save space.
	samtools view -t $reference_genome".fai" -Sb $outputfolder/PoPoolationTE/pe-reads.sam | samtools sort - $outputfolder/PoPoolationTE/pe-reads.sorted
	rm $outputfolder/PoPoolationTE/pe-reads.sam
	samtools view $outputfolder/PoPoolationTE/pe-reads.sorted.bam > $outputfolder/PoPoolationTE/pe-reads.sorted.sam

	# Perform the rest of the PoPoolationTE workflow
	perl identify-te-insertsites.pl --input $outputfolder/PoPoolationTE/pe-reads.sorted.sam --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --narrow-range $readlength --min-count 3 --min-map-qual 15 --output $outputfolder/PoPoolationTE/te-fwd-rev.txt --insert-distance $median_insertsize --read-length $readlength
	perl genomic-N-2gtf.pl --input $reference_genome > $outputfolder/PoPoolationTE/poly_n.gtf
	perl crosslink-te-sites.pl --directional-insertions $outputfolder/PoPoolationTE/te-fwd-rev.txt --min-dist $readlength --max-dist $maxdist --output $outputfolder/PoPoolationTE/te-inserts.txt --single-site-shift 100 --poly-n $outputfolder/PoPoolationTE/poly_n.gtf --te-hierarchy $te_hierarchy --te-hier-level family

	# Convert the gff format to a format that PoPoolationTE understands.
	awk -F'[\t=;]' '{if($7=="-")print $1"\tF\t"$5"\t"$10"\n"$1"\tR\t"$4"\t"$10; else print $1"\tF\t"$4"\t"$10"\n"$1"\tR\t"$5"\t"$10}' $gff_te_locations > $outputfolder/PoPoolationTE/known-te-insertions.txt
	perl update-teinserts-with-knowntes.pl --known $outputfolder/PoPoolationTE/known-te-insertions.txt --output $outputfolder/PoPoolationTE/te-insertions.updated.txt --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --max-dist $maxdist --te-insertions $outputfolder/PoPoolationTE/te-inserts.txt --single-site-shift 100
	perl estimate-polymorphism.pl --sam-file $outputfolder/PoPoolationTE/pe-reads.sorted.sam --te-insert-file $outputfolder/PoPoolationTE/te-insertions.updated.txt --te-hierarchy-file $te_hierarchy --te-hierarchy-level family --min-map-qual 15 --output $outputfolder/PoPoolationTE/te-polymorphism.txt
	perl filter-teinserts.pl --te-insertions $outputfolder/PoPoolationTE/te-polymorphism.txt --output $outputfolder/PoPoolationTE/te-poly-filtered.txt --discard-overlapping --min-count 5

	# Create an output file in BED format selecting only the relevant inserts.
	# Name and description for use with the UCSC genome browser are added to output here.
	echo -e "track name=\"$sample"_PoPoolationTE"\" description=\"$sample"_PoPoolationTE"\"" > $outputfolder/PoPoolationTE/$sample"_popoolationte_raw.bed"
	echo -e "track name=\"$sample"_PoPoolationTE"\" description=\"$sample"_PoPoolationTE"\"" > $outputfolder/PoPoolationTE/$sample"_popoolationte_redundant.bed"
	echo -e "track name=\"$sample"_PoPoolationTE"\" description=\"$sample"_PoPoolationTE"\"" > $outputfolder/PoPoolationTE/$sample"_popoolationte_nonredundant.bed"

	awk -F"\t" -v sample=$sample '{if($9~/-/) {printf "%s\t%d\t%d",$1,$2-1,$16; printf "\t"$13+$20"\t"$4;} else if($16~/-/) {printf "%s\t%d\t%d",$1,$10-1,$2; printf "\t"$13+$20"\t"$4;} else {printf "%s\t%d\t%d",$1,$10-1,$16; printf "\t"$13+$20"\t"$4;} if($7!="-") print "_reference_"sample"_popoolationte_rp_\t0\t.\t"$3"\t"$11"\t"$18; else print "_non-reference_"sample"_popoolationte_rp_\t0\t.\t"$3"\t"$11"\t"$18}' $outputfolder/PoPoolationTE/te-poly-filtered.txt > $outputfolder/PoPoolationTE/$sample"_popoolationte_presort_raw.txt"

	sort -k1,1 -k2,2n $outputfolder/PoPoolationTE/$sample"_popoolationte_presort_raw.txt" > $outputfolder/PoPoolationTE/$sample"_popoolationte_precut_raw.txt"
	awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5NR"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $outputfolder/PoPoolationTE/$sample"_popoolationte_precut_raw.txt" > $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_raw.txt"
	cut -f1-3,5-7 $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_raw.txt" >> $outputfolder/PoPoolationTE/$sample"_popoolationte_raw.bed"

	# Filtering for inserts with support for both ends and more than 10% of reads support an insert for at least one end
	awk '{if($8=="FR" && ($9>=0.1 || $10>=0.1)) print $0}' $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_raw.txt" > $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_redundant.txt"
	cut -f1-3,5-7 $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_redundant.txt" >> $outputfolder/PoPoolationTE/$sample"_popoolationte_redundant.bed"

	# Filtering out redundant insertions
	sort -k1,3 -k4rn $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_redundant.txt" | sort -u -k1,3 | cut -f1-3,5-7 > $outputfolder/PoPoolationTE/tmp
	sort -k1,1 -k2,2n $outputfolder/PoPoolationTE/tmp >> $outputfolder/PoPoolationTE/$sample"_popoolationte_nonredundant.bed"

	# Clean up intermediate files
	rm $outputfolder/PoPoolationTE/tmp $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_raw.txt" $outputfolder/PoPoolationTE/$sample"_popoolationte_counted_precut_redundant.txt" $outputfolder/PoPoolationTE/$sample"_popoolationte_precut_raw.txt" $outputfolder/PoPoolationTE/$sample"_popoolationte_presort_raw.txt" $outputfolder/PoPoolationTE/1.sam $outputfolder/PoPoolationTE/2.sam $outputfolder/PoPoolationTE/pe-reads.sorted.bam

else
	echo "Supply TE masked fasta reference file with consensus TE sequences used to mask added as separate chromosomes as option 1"
	echo "Supply TE hierarchy file as option 2"
	echo "Supply fastq1 as option 3"
	echo "Supply fastq2 as option 4"
	echo "Supply a gff3 annotation of the TEs in the genome as option 5"
	echo "Supply a number of processors as option 6"
	echo "Supply output folder as option 7"
fi

