#!/bin/bash

usage () 
{
    echo "-s : sample name [Required]"
    echo "-R : referencefolder [Required]"
    echo "-o : An output folder for the run. [Optional]"
    echo "-1 : The absolute path of the first fastq file from a paired end read, this should be named ending _1.fastq. [Required]"
	echo "-2 : The absolute path of the second fastq file from a paired end read, this should be named ending _2.fastq. [Optional]"
	echo "-r : A reference genome sequence in fasta format. [Required]"
	echo "-c : The consensus sequences of the TEs for the species in fasta format. [Required]"
    echo "-m : McC location [Required]"
    echo "-p : The number of processors to use for parallel stages of the pipeline. [Optional: default = 1]"
	echo "-i : If this option is specified then all sample specific intermediate files will be removed, leaving only"
	echo "     the overall results. The default is to leave sample specific intermediate files (may require large amounts"
	echo "     of disk space)"
	echo "-h : Prints this help guide."
}

# Set default value for processors in case it is not supplied
processors=1

# Get the options supplied to the program
while getopts ":s:R:o:1:2:r:c:m:p:hi" opt;
do
	case $opt in
		s)
			sample=$OPTARG
			;;
		R)
			referencefolder=$OPTARG
			;;
		o)
			te_cov_dir=$OPTARG
			;;
		1)
			fastq1=$OPTARG
			;;
		2)	
			fastq2=$OPTARG
			;;
		r)
			reference_genome=$OPTARG
			;;
		c)
			consensus_te_seqs=$OPTARG
			;;
        m)
            mcclintock_location=$OPTARG
            ;;
        p)
            processors=$OPTARG
            ;;
        i)
            remove_intermediates=on
            ;;
		h)
			usage
			exit 1
			;;
		\?)
			echo "Unknown option: -$OPTARG"
			usage
			exit 1
			;;
		:)
			echo "Missing option argument for -$OPTARG"
			usage
			exit 1
			;;
	esac
done

# create tmp folder for intermediate files
tmp_dir=$te_cov_dir/te_cov_tmp
mkdir -p $tmp_dir
rm -rf $tmp_dir/*

genome=${reference_genome##*/}
genome=${reference_genome%%.*}

# Hard mask reference genome to exclude reference TE sequences using repeatmasker
if [[ ! -f $reference_genome".masked" || ! -f $reference_genome".out.gff" ]]
then
    time RepeatMasker -dir $referencefolder -s -gff -nolow -no_is -a $reference_genome -lib $consensus_te_seqs -pa $processors
    # RepeatMasker appears to override the custom database names during the ProcessRepeats step so this changes them back for
    # Drosophila, more rules like this may be needed for other reference genomes
    sed "s/McClintock-int/McClintock/g" $reference_genome".out.gff" > $referencefolder/tmp
    sed "s/POGON1/pogo/g" $referencefolder/tmp > $reference_genome".out.gff"
    perl $mcclintock_location/scripts/fixfastalinelength.pl $reference_genome".masked" 80 $reference_genome".masked2"
    mv $reference_genome".masked2" $reference_genome".masked"
fi
ref_masked=$reference_genome".masked"

# create augmented genome by including TE sequences at the end of the masked reference genome.
if [[ ! -f $ref_masked".aug" ]]
then
	cat $ref_masked $consensus_te_seqs > $ref_masked".aug"
fi
ref_masked_aug=$ref_masked".aug"

# Create indexes of augmented masked reference genome if not already made
if [[ ! -f $ref_masked_aug".fai" ]]
then
	samtools faidx $ref_masked_aug
fi
if [[ ! -f $ref_masked_aug".bwt" ]]
then
	bwa index $ref_masked_aug
fi

# Map reads to the augmented genome.
bwa mem -t $processors -v 0 -R '@RG\tID:'$sample'\tSM:'$sample $ref_masked_aug $fastq1 $fastq2 > $tmp_dir/$sample.sam
sam=$tmp_dir/$sample.sam
samtools view -Sb -t $ref_masked_aug".fai" $sam | samtools sort - $tmp_dir/$sample
bam=$tmp_dir/$sample.bam
samtools index $bam

# get non-TE regions of the reference genome in bed format
if [[ ! -f $referencefolder/$genome".fasta.out.complement.bed" ]]
then
	rep_out=$reference_genome".out"
	rep_bed=$reference_genome".bed"
	rep_bed_sorted=$reference_genome".sorted.bed"
	awk 'BEGIN{OFS="\t"}{if(NR>3) {if($9=="C"){strand="-"}else{strand="+"};print $5,$6-1,$7,$10,".",strand}}' $rep_out > $rep_bed
	bedtools sort -i $rep_bed > $rep_bed_sorted
	grep chr $ref_masked_aug.fai | cut -f1,2 | sort > $referencefolder/$genome".sorted.genome"
	bedtools slop -i $rep_bed_sorted -g $referencefolder/$genome".sorted.genome" -l 1 -r 0 > $referencefolder/$genome".fasta.out.zero.bed"
	bedtools complement -i $referencefolder/$genome".fasta.out.zero.bed" -g $referencefolder/$genome".sorted.genome" > $referencefolder/$genome".fasta.out.complement.bed"
fi
bed_nonte=$referencefolder/$genome".fasta.out.complement.bed"

# get TE list from consensus sequence file
if [[ ! -f $referencefolder/te_list ]]
then
	grep ">" $consensus_te_seqs | sed 's/>//g' > $referencefolder/te_list
fi
te_list=$referencefolder/te_list

# Calculate depth of coverage for every TE using Samtools depth module on the BAM file and normalize by average coverage in no TE regione of the normal reference genome.
genome_avg_depth=`samtools depth -b  $bed_nonte $bam | awk '{ total += $3 } END { print total/NR }'`
echo $genome_avg_depth > $referencefolder/genome_avg_depth
printf '%s\n' "TE Family" "Normalized Depth" | paste -sd ',' > $te_cov_dir/te_depth.csv
for te in `cat $te_list`
do
    te_depth=`samtools depth -r $te $bam | awk '{ total += $3 } END { print total/NR }'`
    te_depth_normalized=$(echo "$te_depth / $genome_avg_depth" | bc -l )
    printf '%s,%.2f\n' "$te" "$te_depth_normalized" | paste -sd ',' >> $te_cov_dir/te_depth.csv
done

# # remove tmp folder
# if [[ "$remove_intermediates" = "on" ]]
# then
#     rm -rf $tmp_dir
# fi