#! /bin/bash -l

usage () 
{
	echo "McClintock Usage"
	echo "This script takes the following inputs and will run 5 different transposable element (TE) detection methods:"
	echo "-r : A reference genome sequence in fasta format. [Required]"
	echo "-c : The consensus sequences of the TEs for the species in fasta format. [Required]"
	echo "-g : The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID"
	echo "     attribute for every entry. [Required]"
	echo "-t : A tab delimited file with one entry per ID in the GFF file and two columns: the first containing"
	echo "     the ID and the second containing the TE family it belongs to. The family should correspond to the"
	echo "     names of the sequences in the consensus fasta file. [Required]"
	echo "-b : Retain the sorted and indexed BAM file of the paired end data aligned to the reference genome."
	echo "-i : If this option is specified then all sample specific intermediate files will be removed, leaving only"
	echo "     the overall results. The default is to leave sample specific intermediate files (may require large amounts"
	echo "     of disk space)"
	echo "-1 : The absolute path of the first fastq file from a paired end read, this should be named ending _1.fastq. [Required]"
	echo "-2 : The absolute path of the second fastq file from a paired end read, this should be named ending _2.fastq. [Required]"
	echo "-p : The number of processors to use for parallel stages of the pipeline."
	echo "-h : Prints this help guide."
}

# Set default value for processors in case it is not supplied
processors=1

# Get the options supplied to the program
while getopts ":r:c:g:t:1:2:p:hib" opt;
do
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
		p)
			processors=$OPTARG
			;;
		i)	
			remove_intermediates=on
			;;
		b)
			save_bam=on
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

# Test for presence of required arguments
if [[ -z "$inputr" || -z "$inputc" || -z "$inputg" || -z "$inputt" || -z "$input1" || -z "$input2" ]];
then
	echo "A required parameter is missing"
	usage
	exit 1
fi

# Set up folder structure

printf "\nCreating directory structure...\n\n" | tee -a /dev/stderr

genome=${inputr##*/}
genome=${genome%%.*}
sample=${input1##*/}
sample=${sample%%_1.f*}

test_dir=`pwd`
if [ ! -d $test_dir/$genome ];
then
	mkdir $test_dir/$genome/
	mkdir $test_dir/$genome/reference
fi
mkdir $test_dir/$genome/$sample
mkdir $test_dir/$genome/$sample/reads
mkdir $test_dir/$genome/$sample/bam
mkdir $test_dir/$genome/$sample/sam
mkdir $test_dir/$genome/$sample/results
mkdir $test_dir/$genome/$sample/results/qualitycontrol
mkdir $test_dir/$genome/$sample/results/originalmethodresults

# Copy input files in to sample directory (neccessary for RelocaTE)
reference_genome_file=${inputr##*/}
if [ ! -f $test_dir/$genome/reference/$reference_genome_file ];
then
	cp -n $inputr $test_dir/$genome/reference/$reference_genome_file
fi
consensus_te_seqs_file=${inputc##*/}
if [ ! -f $test_dir/$genome/reference/$consensus_te_seqs_file ];
then
	cp -n $inputc $test_dir/$genome/reference/$consensus_te_seqs_file
fi
te_locations_file=${inputg##*/}
if [ ! -f $test_dir/$genome/reference/$te_locations_file ];
then
	cp -n $inputg $test_dir/$genome/reference/$te_locations_file
fi
te_families_file=${inputt##*/}
if [ ! -f $test_dir/$genome/reference/$te_families_file ];
then
	cp -n $inputt $test_dir/$genome/reference/$te_families_file
fi
fastq1_file=${input1##*/}
cp -s $input1 $test_dir/$genome/$sample/reads/$fastq1_file
fastq2_file=${input2##*/}
cp -s $input2 $test_dir/$genome/$sample/reads/$fastq2_file

# Assign variables to input files
reference_genome=$test_dir/$genome/reference/$reference_genome_file
consensus_te_seqs=$test_dir/$genome/reference/$consensus_te_seqs_file
te_locations=$test_dir/$genome/reference/$te_locations_file
te_families=$test_dir/$genome/reference/$te_families_file
fastq1=$test_dir/$genome/$sample/reads/$fastq1_file
fastq2=$test_dir/$genome/$sample/reads/$fastq2_file

# If fastqc is installed then launch fastqc on the input fastqs
location=`which fastqc`
if test -z "$location"
then
	printf "\nfastqc not installed, skipping input quality analysis...\n\n" | tee -a /dev/stderr
else
	mkdir $test_dir/$genome/$sample/results/qualitycontrol/fastqc_analysis
	fastqc -t $processors $fastq1 $fastq2 -o $test_dir/$genome/$sample/results/qualitycontrol/fastqc_analysis
fi

# Create indexes of reference genome if not already made for this genome
if [ ! -f $reference_genome".fai" ];
then
	samtools faidx $reference_genome
fi
if [ ! -f $reference_genome".bwt" ];
then
	bwa index $reference_genome
fi

# Extract sequence of all reference TE copies if this has not already been done
# Cut first line if it begins with #
if [ ! -f $test_dir/$genome/reference/all_te_seqs.fasta ];
then
	grep -v '^#' $te_locations | awk -F'[\t=;]' 'BEGIN {OFS = "\t"}; {printf $1"\t"$2"\t"; for(x=1;x<=NF;x++) if ($x~"ID") printf $(x+1); print "\t"$4,$5,$6,$7,$8,"ID="}' | awk -F'\t' '{print $0$3";Name="$3";Alias="$3}' > edited.gff
	mv edited.gff $te_locations
	bedtools getfasta -name -fi $reference_genome -bed $te_locations -fo $test_dir/$genome/reference/all_te_seqs.fasta
fi
all_te_seqs=$test_dir/$genome/reference/all_te_seqs.fasta

# Create bed file of reference TE locations
if [ ! -f $test_dir"/"$genome"/reference/reference_TE_locations.bed" ];
then
	awk -F["\t"\;=] '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' $te_locations > $test_dir/$genome/reference/reference_TE_locations.bed
fi
bed_te_locations_file=$test_dir/$genome/reference/reference_TE_locations.bed

# Create sam files for input
printf "\nCreating bam alignment...\n\n" | tee -a /dev/stderr

bwa mem -a -t $processors -v 0 $reference_genome $fastq1 $fastq2 > $test_dir/$genome/$sample/sam/$sample.sam
sort --temporary-directory=. $test_dir/$genome/$sample/sam/$sample.sam > $test_dir/$genome/$sample/sam/sorted$sample.sam
rm $test_dir/$genome/$sample/sam/$sample.sam
mv $test_dir/$genome/$sample/sam/sorted$sample.sam $test_dir/$genome/$sample/sam/$sample.sam 
sam=$test_dir/$genome/$sample/sam/$sample.sam
sam_folder=$test_dir/$genome/$sample/sam

# Calculate the median insert size of the sample
median_insertsize=`cut -f9 $sam | sort -n | awk '{if ($1 > 0) ins[reads++]=$1; } END { print ins[int(reads/2)]; }'`
printf "\nMedian insert size = $median_insertsize\n\n" | tee -a /dev/stderr
echo $median_insertsize > $test_dir/$genome/$sample/results/qualitycontrol/median_insertsize

# Create bam files for input
samtools view -Sb $sam > $test_dir/$genome/$sample/bam/$sample.bam
samtools sort $test_dir/$genome/$sample/bam/$sample.bam $test_dir/$genome/$sample/bam/sorted$sample
rm $test_dir/$genome/$sample/bam/$sample.bam
mv $test_dir/$genome/$sample/bam/sorted$sample.bam $test_dir/$genome/$sample/bam/$sample.bam
bam=$test_dir/$genome/$sample/bam/$sample.bam
samtools index $bam

# Get stats of bam file from samtools
samtools flagstat $bam > $test_dir/$genome/$sample/results/qualitycontrol/bwamem_bamstats.txt

################################## Run TE-locate ##################################

printf "\nRunning TE-locate pipeline...\n\n" | tee -a /dev/stderr

# Adjust hierachy levels
cd TE-locate
telocate_te_locations=${te_locations%.*}
telocate_te_locations=$telocate_te_locations"_HL.gff"
if [ ! -f $telocate_te_locations ];
then
	perl TE_hierarchy.pl $te_locations $te_families Alias
fi

bash runtelocate.sh $sam_folder $reference_genome $telocate_te_locations 2 $sample $median_insertsize

# Save the original result file and the bed files filtered by mcclintock
mv $sample/$sample"_telocate"* $test_dir/$genome/$sample/results/
mkdir $test_dir/$genome/$sample/results/originalmethodresults/TE-locate
cp $sample/*.info $test_dir/$genome/$sample/results/originalmethodresults/TE-locate

# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [ "$remove_intermediates" = "on" ]
then
	printf "\nRemoving TE-locate intermediate files\n\n"
	rm -r $sam_folder
	rm -r $sample
fi

################################## Run RetroSeq ##################################

printf "\nRunning RetroSeq pipeline...\n\n" | tee -a /dev/stderr

cd ../RetroSeq
bash runretroseq.sh $consensus_te_seqs $bam $reference_genome $bed_te_locations_file $te_families

# Save the original result file and the bed files filtered by mcclintock
mv $sample/$sample"_retroseq"* $test_dir/$genome/$sample/results/
mkdir $test_dir/$genome/$sample/results/originalmethodresults/RetroSeq
cp $sample/$sample".calling.PE.vcf" $test_dir/$genome/$sample/results/originalmethodresults/RetroSeq

# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [ "$remove_intermediates" = "on" ]
then
	printf "\nRemoving RetroSeq intermediate files\n\n"
	# If the save bam option is specified then override the command to delete the bam files.
	if [ "$save_bam" != "on" ]
	then
		rm -r $test_dir/$genome/$sample/bam
	fi
	rm -r $sample
fi

################################## Run TEMP ##################################

printf "\nRunning TEMP pipeline...\n\n" | tee -a /dev/stderr

cd ../TEMP
bash runtemp.sh $reference_genome $fastq1 $fastq2 $consensus_te_seqs $bed_te_locations_file $te_families $median_insertsize $sample $processors

# Save the original result file and the bed files filtered by mcclintock
mv $sample/$sample"_temp"* $test_dir/$genome/$sample/results/
mkdir $test_dir/$genome/$sample/results/originalmethodresults/TEMP
cp $sample/$sample".insertion.refined.bp.summary" $test_dir/$genome/$sample/results/originalmethodresults/TEMP

# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [ "$remove_intermediates" = "on" ]
then
	printf "\nRemoving TEMP intermediate files\n\n"
	rm -r $sample
fi

################################## Run RelocaTE ##################################

printf "\nRunning RelocaTE pipeline...\n\n" | tee -a /dev/stderr

# Add TSD lengths to consensus TE sequences
if [ ! -f $test_dir/$genome/reference/relocate_te_seqs.fasta ];
then
	awk '{if (/>/) print $0" TSD=UNK"; else print $0}' $consensus_te_seqs > $test_dir/$genome/reference/relocate_te_seqs.fasta
fi
relocate_te_seqs=$test_dir/$genome/reference/relocate_te_seqs.fasta

# Create general gff file to allow reference TE detection in RelocaTE
if [ ! -f $test_dir/$genome/reference/relocate_te_locations.gff ];
then
	awk 'FNR==NR{array[$1]=$2;next}{print $1,$2,array[$3],$4,$5,$6,$7,$8,$9}' FS='\t' OFS='\t' $te_families $te_locations > $test_dir/$genome/reference/relocate_te_locations.gff
fi
relocate_te_locations=$test_dir/$genome/reference/relocate_te_locations.gff

cd ../RelocaTE
bash runrelocate.sh $relocate_te_seqs $reference_genome $test_dir/$genome/$sample/reads $sample $relocate_te_locations

# Save the original result file and the bed files filtered by mcclintock
mv $sample/$sample"_relocate"* $test_dir/$genome/$sample/results/
mkdir $test_dir/$genome/$sample/results/originalmethodresults/RelocaTE
cp -r $sample/*/results $test_dir/$genome/$sample/results/originalmethodresults/RelocaTE

# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [ "$remove_intermediates" = "on" ]
then
	printf "\nRemoving RelocaTE intermediate files\n\n"
	rm -r $sample
fi

################################## Run ngs_te_mapper ##################################

printf "\nRunning ngs_te_mapper pipeline...\n\n" | tee -a /dev/stderr

cd ../ngs_te_mapper

bash runngstemapper.sh $consensus_te_seqs $reference_genome $sample $fastq1 $fastq2 

# Save the original result file and the bed file filtered by mcclintock
mv $sample/$sample"_ngs_te_mapper_nonredundant.bed" $test_dir/$genome/$sample/results/
mkdir $test_dir/$genome/$sample/results/originalmethodresults/ngs_te_mapper
cp $sample/analysis/bed_tsd/*.bed $test_dir/$genome/$sample/results/originalmethodresults/ngs_te_mapper

# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [ "$remove_intermediates" = "on" ]
then
	printf "\nRemoving ngs_te_mapper intermediate files\n\n"
	rm -r $sample
fi

################################## Run PoPoolationTE ##################################

printf "\nRunning PoPoolationTE pipeline...\n\n" | tee -a /dev/stderr

# Create te_hierachy
if [ ! -f $test_dir/$genome/reference/te_hierarchy ]; then
	printf "insert\tid\tfamily\tsuperfamily\tsuborder\torder\tclass\tproblem\n" > $test_dir/$genome/reference/te_hierarchy
	awk '{printf $0"\t"$2"\t"$2"\tna\tna\tna\t0\n"}' $te_families >> $test_dir/$genome/reference/te_hierarchy
fi
te_hierarchy=$test_dir/$genome/reference/te_hierarchy

cd ../PoPoolationTE
bash runpopoolationte.sh $reference_genome $all_te_seqs $te_hierarchy $fastq1 $fastq2 $te_locations $processors

# Save the original result file and the bed files filtered by mcclintock
mv $sample/$sample"_popoolationte"* $test_dir/$genome/$sample/results/
mkdir $test_dir/$genome/$sample/results/originalmethodresults/PoPoolationTE
cp $sample/te-poly-filtered.txt $test_dir/$genome/$sample/results/originalmethodresults/PoPoolationTE

# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [ "$remove_intermediates" = "on" ]
then
	printf "\nRemoving PoPoolationTE intermediate files\n\n"
	rm -r $sample
fi

#########################################################################################

cd ..
# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [ "$remove_intermediates" = "on" ]
then
	printf "\nRemoving McClintock intermediate files\n\n"
	rm -r $genome/$sample/reads
fi

printf "\nPipeline Complete\n\n" | tee -a /dev/stderr
