#!/bin/bash -l

# Get the options supplied to the program.
while getopts ":r:c:g:t:1:2:h" opt; do
	case $opt in
		r)
			inputr=$OPTARG
			;;
		c)
			inputc=$OPTARG
			;;
		g)
			inputg=$OPTARG
			;;
		t)
			inputt=$OPTARG
			;;
		1)	
			input1=$OPTARG
			;;
		2)
			input2=$OPTARG
			;;
		h)
			echo "This script takes the following inputs and will run 5 different transposable element (TE) detection methods:"
			echo "-r : A reference genome sequence in fasta format."
			echo "-c : The consensus sequences of the TEs for the species in fasta format."
			echo "-g : The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry."
			echo "-t : A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file."
			echo "-1 : The absolute path of the first fastq file from a paired end read, this should be named ending _1.fastq."
			echo "-2 : The absolute path of the second fastq file from a paired end read, this should be named ending _2.fastq."
			echo "-h : Prints this help guide."
			exit 1
			;;
		\?)
			echo "Invalid option: -$OPTARG"
			echo "This script takes the following inputs and will run 5 different transposable element (TE) detection methods:"
			echo "-r : A reference genome sequence in fasta format."
			echo "-c : The consensus sequences of the TEs for the species in fasta format."
			echo "-g : The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry."
			echo "-t : A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file."
			echo "-1 : The absolute path of the first fastq file from a paired end read, this should be named ending _1.fastq."
			echo "-2 : The absolute path of the second fastq file from a paired end read, this should be named ending _2.fastq."
			echo "-h : Prints this help guide."
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument."
			echo "This script takes the following inputs and will run 5 different transposable element (TE) detection methods:"
			echo "-r : A reference genome sequence in fasta format."
			echo "-c : The consensus sequences of the TEs for the species in fasta format."
			echo "-g : The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry."
			echo "-t : A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file."
			echo "-1 : The absolute path of the first fastq file from a paired end read, this should be named ending _1.fastq."
			echo "-2 : The absolute path of the second fastq file from a paired end read, this should be named ending _2.fastq."
			echo "-h : Prints this help guide.\n"
			exit 1
			;;
	esac
done

# Set up folder structure

printf "\nCreating directory structure...\n\n"

genome=${inputr##*/}
genome=${genome%%.*}
sample=${input1##*/}
sample=${sample%%_1.f*}

test_dir=`pwd`
mkdir $test_dir/$genome/
mkdir $test_dir/$genome/reference
mkdir $test_dir/$genome/$sample
mkdir $test_dir/$genome/$sample/fastq
mkdir $test_dir/$genome/$sample/bam
mkdir $test_dir/$genome/$sample/sam

# Copy inout files in to sample directory (neccessary for RelocaTE)
reference_genome_file=${inputr##*/}
cp -n $inputr $test_dir/$genome/reference/$reference_genome_file
consensus_te_seqs_file=${inputc##*/}
cp -n $inputc $test_dir/$genome/reference/$consensus_te_seqs_file
te_locations_file=${inputg##*/}
cp -n $inputg $test_dir/$genome/reference/$te_locations_file
te_families_file=${inputt##*/}
cp -n $inputt $test_dir/$genome/reference/$te_families_file
fastq1_file=${input1##*/}
cp -s $input1 $test_dir/$genome/$sample/fastq/$fastq1_file
fastq2_file=${input2##*/}
cp -s $input2 $test_dir/$genome/$sample/fastq/$fastq2_file

# Assign variables to input files
reference_genome=$test_dir/$genome/reference/$reference_genome_file
consensus_te_seqs=$test_dir/$genome/reference/$consensus_te_seqs_file
te_locations=$test_dir/$genome/reference/$te_locations_file
te_families=$test_dir/$genome/reference/$te_families_file
fastq1=$test_dir/$genome/$sample/fastq/$fastq1_file
fastq2=$test_dir/$genome/$sample/fastq/$fastq2_file

# Create indexes of reference genome
samtools faidx $reference_genome
bwa index $reference_genome

# Extract sequence of all reference TE copies
# Cut first line if it begins with #
grep -v '^#' $te_locations | awk -F'[\t=;]' 'BEGIN {OFS = "\t"}; {printf $1"\t"$2"\t"; for(x=1;x<=NF;x++) if ($x~"ID") printf $(x+1); print "\t"$4,$5,$6,$7,$8,"ID="}' | awk -F'\t' '{print $0$3";Name="$3";Alias="$3}' > edited.gff
mv edited.gff $te_locations
bedtools getfasta -name -fi $reference_genome -bed $te_locations -fo $test_dir/$genome/reference/all_te_seqs.fasta
all_te_seqs=$test_dir/$genome/reference/all_te_seqs.fasta

# Create sam and bam files for input

printf "\nCreating bam alignment...\n\n"

bwa mem -v 0 $reference_genome $fastq1 $fastq2 > $test_dir/$genome/$sample/sam/$sample.sam
sort --temporary-directory=. $test_dir/$genome/$sample/sam/$sample.sam > $test_dir/$genome/$sample/sam/sorted$sample.sam
rm $test_dir/$genome/$sample/sam/$sample.sam
mv $test_dir/$genome/$sample/sam/sorted$sample.sam $test_dir/$genome/$sample/sam/$sample.sam 
sam=$test_dir/$genome/$sample/sam/$sample.sam
sam_folder=$test_dir/$genome/$sample/sam

samtools view -Sb $sam > $test_dir/$genome/$sample/bam/$sample.bam
samtools sort $test_dir/$genome/$sample/bam/$sample.bam $test_dir/$genome/$sample/bam/sorted$sample
rm $test_dir/$genome/$sample/bam/$sample.bam 
mv $test_dir/$genome/$sample/bam/sorted$sample.bam $test_dir/$genome/$sample/bam/$sample.bam 
bam=$test_dir/$genome/$sample/bam/$sample.bam 
samtools index $bam

# Run RelocaTE

printf "\nRunning RelocaTE pipeline...\n\n"

# Add TSD lengths to consensus TE sequences
awk '{if (/>/) print $0" TSD=UNK"; else print $0}' $consensus_te_seqs > $test_dir/$genome/reference/relocate_te_seqs.fasta
relocate_te_seqs=$test_dir/$genome/reference/relocate_te_seqs.fasta

# Create general gff file to allow reference TE detection in RelocaTE
awk 'FNR==NR{array[$1]=$2;next}{print $1,$2,array[$3],$4,$5,$6,$7,$8,$9}' FS='\t' OFS='\t' $te_families $te_locations > $test_dir/$genome/reference/relocate_te_locations.gff
relocate_te_locations=$test_dir/$genome/reference/relocate_te_locations.gff

cd RelocaTE
bash runrelocate.sh $relocate_te_seqs $reference_genome $test_dir/$genome/$sample/fastq $sample $relocate_te_locations

# Run ngs_te_mapper pipeline

printf "\nRunning ngs_te_mapper pipeline...\n\n"

cd ../ngs_te_mapper

bash runngstemapper.sh $consensus_te_seqs $reference_genome $sample $fastq1 $fastq2 

# Run RetroSeq

printf "\nRunning RetroSeq pipeline...\n\n"

cd ../RetroSeq
bash runretroseq.sh $consensus_te_seqs $bam $reference_genome 

# Run TE-locate

printf "\nRunning TE-locate pipeline...\n\n"

# Adjust hierachy levels
cd ../TE-locate
perl TE_hierarchy.pl $te_locations $te_families Alias
telocate_te_locations=${te_locations%.*}
telocate_te_locations=$telocate_te_locations"_HL.gff"

bash runtelocate.sh $sam_folder $reference_genome $telocate_te_locations 2 $sample

# Run PoPoolationTE

printf "\nRunning PoPoolationTE pipeline...\n\n"

# Create te_hierachy
printf "insert\tid\tfamily\tsuperfamily\tsuborder\torder\tclass\tproblem\n" > $test_dir/$genome/reference/te_hierarchy
awk '{printf $0"\t"$2"\t"$2"\tna\tna\tna\t0\n"}' $te_families >> $test_dir/$genome/reference/te_hierarchy
te_hierarchy=$test_dir/$genome/reference/te_hierarchy

cd ../popoolationte
bash runpopoolationte.sh $reference_genome $all_te_seqs $te_hierarchy $fastq1 $fastq2 $te_locations
