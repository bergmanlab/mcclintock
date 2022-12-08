#! /bin/bash -l

usage () 
{
	printf "McClintock version "
	git rev-parse HEAD
	echo "This script takes the following inputs and will run 5 different transposable element (TE) detection methods:"
	echo "-r : A reference genome sequence in fasta format. [Required]"
	echo "-c : The consensus sequences of the TEs for the species in fasta format. [Required]"
	echo "-g : The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID"
	echo "     attribute for every entry. [Optional]"
	echo "-t : A tab delimited file with one entry per ID in the GFF file and two columns: the first containing"
	echo "     the ID and the second containing the TE family it belongs to. The family should correspond to the"
	echo "     names of the sequences in the consensus fasta file. [Optional - required if GFF (option -g) is supplied]"
	echo "-1 : The absolute path of the first fastq file from a paired end read, this should be named ending _1.fastq. [Required]"
	echo "-2 : The absolute path of the second fastq file from a paired end read, this should be named ending _2.fastq. [Optional]"
	echo "-o : An output folder for the run. This should be the absolute path. If not supplied then the directory mcclintock"
	echo "     is launched from will be used. [Optional]"
	echo "-T : If this option is specified then fastq comments (e.g. barcode) will be incorporated to SAM output. Warning: do not"
	echo "     use this option if the input fastq files do not have comments. [Optional: default will not include this]"
	echo "-d : If this option is specified then McClintock will perform depth of coverage analysis for every TE. Note: Doing"
	echo "     TE-based coverage analysis will result in longer running time. A fasta file can be provided here for coverage"
	echo "     analysis. If no file is provided here, the consensus sequences of the TEs will be used for the analysis. [Optional:"
	echo "     default will not include this]"
	echo "-D : IF this option is specified then only depth of coverage analysis for TEs will be performed. [Optional]"
	echo "-s : A fasta file that will be used for TE-based coverage analysis, if not supplied then the consensus sequences"
	echo "     of the TEs will be used for the analysis. [Optional]"
	echo "-b : Retain the sorted and indexed BAM file of the paired end data aligned to the reference genome."
	echo "-i : If this option is specified then all sample specific intermediate files will be removed, leaving only"
	echo "     the overall results. The default is to leave sample specific intermediate files (may require large amounts"
	echo "     of disk space)"
	echo "-C : This option will include the consensus TE sequences as extra chromosomes in the reference file (useful if the "
	echo "     organism is known to have TEs that are not present in the reference strain). [Optional: default will not include"
	echo "     this]"
	echo "-R : This option will include the reference TE sequences as extra chromosomes in the reference file [Optional: default"
	echo "     will not include this]"
	echo "-m : A string containing the list of software you want the pipeline to use for analysis e.g. \"-m relocate TEMP "
	echo "     ngs_te_mapper\" will launch only those three methods [Optional: default is to run all methods]"
	echo "-p : The number of processors to use for parallel stages of the pipeline. [Optional: default = 1]"
	echo "-M : The amount of memory available for the pipeline in GB. [Optional: default = 4]"
	echo "-h : Prints this help guide."
}

# set up conda env for base mcclintock pipeline
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate MCCLINTOCK_main

# Set default value for processors in case it is not supplied
processors=1
memory=4
# Default behaviour is to run all methods if no option is supplied
methods="ngs_te_mapper RelocaTE TEMP RetroSeq PoPoolationTE TE-locate"
# If an output folder is not specified then default to th current directory
outputfolder=`pwd`
# Save the location of mcclintock
mcclintock_location="$( cd "$(dirname "$0")" ; pwd -P )"
cd $mcclintock_location/

# Get the options supplied to the program
while getopts ":r:c:g:t:s:1:2:o:p:M:m:hiTdDbCR" opt;
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
		s)
			inputff=$OPTARG
			;;
		1)	
			input1=$OPTARG
			;;
		2)
			input2=$OPTARG
			;;
		o)
			outputfolder=$OPTARG
			;;
		p)
			processors=$OPTARG
			;;
		M)
			memory=$OPTARG
			;;
		m)
			methods=$OPTARG
			;;
		C)
			addconsensus=on
			;;
		R)
			addrefcopies=on
			;;
		i)
			remove_intermediates=on
			;;
		T)
			save_tag=on
			;;
		d)
			te_cov=on
			;;
		D)
			cov_only=on
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
if [[ -z "$inputr" || -z "$inputc" || -z "$input1" ]]
then
	echo "A required parameter is missing"
	usage
	exit 1
fi

# If a GFF file is supplied then a TE family file that links it to the fasta consensus is also needed
if [[ "$inputg" ]]
then
	if [[ -z "$inputt" ]]
	then
		echo "If a GFF file is supplied then a TE family file that links it to the fasta consensus is also needed"
		usage
		exit 1
	fi
fi

# Output Version number and date

date=`date +%d_%m_%y`
printf "\nRunning McClintock version: " | tee -a /dev/stderr
git rev-parse HEAD | tee -a /dev/stderr
printf "\n\nDate of run is $date\n\n" | tee -a /dev/stderr

# If only one fastq has been supplied assume this is single ended data and launch only ngs_te_mapper and RelocaTE
single_end="false"
if [[ "$input1" && -z "$input2" ]]
then
	echo "Assuming single ended data and launching only ngs_te_mapper and RelocaTE"
	methods="ngs_te_mapper RelocaTE"
	single_end="true"
fi

# If coverage only option is provided then launch only coverage module
if [[ "$cov_only" == "on" ]]
then
	echo "Coverage only option is provided and launching only coverage module"
	methods=""
fi

genome=${inputr##*/}
genome=${genome%%.*}

if [[ $single_end != "true" ]]
then
	sample=${input1##*/}
	sample=${sample%%_1.f*}
	sample=${sample%%.f*}
else
	sample=${input1##*/}
	sample=${sample%%.f*}
fi

# Set up folder structure
printf "\nCreating directory structure...\n\n" | tee -a /dev/stderr

if [[ ! -d $outputfolder/$genome ]]
then
	mkdir -p $outputfolder/
	mkdir $outputfolder/$genome/
	mkdir $outputfolder/$genome/reference
fi
samplefolder=$outputfolder/$genome/$sample
referencefolder=$outputfolder/$genome/reference
mkdir $samplefolder
mkdir $samplefolder/reads
mkdir $samplefolder/bam
mkdir $samplefolder/sam
mkdir $samplefolder/results
mkdir $samplefolder/results/summary
mkdir $samplefolder/results/originalmethodresults

# Copy input files in to sample directory
# Copy the reference fasta file to the run folder
reference_genome_file=${inputr##*/}
if [[ ! -f $referencefolder/$reference_genome_file ]]
then
	# Use script to fix the line length of reference input to 80 characters (needed for samtools index)
	perl $mcclintock_location/scripts/fixfastalinelength.pl $inputr 80 $referencefolder/$reference_genome_file
	samtools faidx $referencefolder/$reference_genome_file
	grep \> $referencefolder/$reference_genome_file | cut -d\> -f2 > $referencefolder/chromosome_names
fi
reference_genome=$referencefolder/$reference_genome_file
# normal_ref_genome=$referencefolder/$reference_genome_file.old
# cp $reference_genome $normal_ref_genome
# samtools faidx $normal_ref_genome
chromosome_names=$referencefolder/chromosome_names

# Copy the TE consesnus fasta file to the run folder
consensus_te_seqs_file=${inputc##*/}
if [[ ! -f $referencefolder/$consensus_te_seqs_file ]]
then
	# Use script to fix the line length of reference input to 80 characters (needed for samtools index)
	perl $mcclintock_location/scripts/fixfastalinelength.pl $inputc 80 $referencefolder/$consensus_te_seqs_file
fi
consensus_te_seqs=$referencefolder/$consensus_te_seqs_file

# Create symbolic links for fastq files to save time and space
fastq1_file=${input1##*/}
cp -s $input1 $samplefolder/reads/$fastq1_file
fastq1=$samplefolder/reads/$fastq1_file
if [[ $single_end != "true" ]]
then
	fastq2_file=${input2##*/}
	cp -s $input2 $samplefolder/reads/$fastq2_file
	fastq2=$samplefolder/reads/$fastq2_file
else
	cp -s $input1 $samplefolder/reads/$sample.unPaired.fastq
fi

# Coverage analysis for each TE
if [[ "$te_cov" == "on" || "$cov_only" == "on" ]]
then
	printf "\nCalculating normalized average coverage for TEs...\n\n" | tee -a /dev/stderr
	if [[ "$inputff" ]]
	then
		printf "\nCustom TE library will be used for coverage analysis.\n\n" | tee -a /dev/stderr
		cov_fasta_file=${inputff##*/}
		cov_fasta_file=${cov_fasta_file%%.*}_cov.fasta
		if [[ ! -f $referencefolder/$cov_fasta_file ]]
		then
			perl $mcclintock_location/scripts/fixfastalinelength.pl $inputff 80 $referencefolder/$cov_fasta_file
		fi
		cov_fasta=$referencefolder/$cov_fasta_file
	else
		printf "\nThe consensus sequences of the TEs will be used for coverage analysis.\n\n" | tee -a /dev/stderr
		cov_fasta=$consensus_te_seqs
	fi
	te_cov_dir=$samplefolder/results/te_coverage
	mkdir -p $te_cov_dir
	ref_genome=$referencefolder/$reference_genome_file

	# load conda env for coverage analysis
	conda deactivate
	conda activate MCCLINTOCK_cov
	if [[ $single_end == "true" ]]
	then
		if [[ "$remove_intermediates" == "on" ]]
		then
			bash $mcclintock_location/scripts/te_coverage.sh -s $sample -R $referencefolder -o $te_cov_dir -1 $fastq1 -r $ref_genome -c $cov_fasta -m $mcclintock_location -p $processors -i -C $chromosome_names -g $genome
		else
			bash $mcclintock_location/scripts/te_coverage.sh -s $sample -R $referencefolder -o $te_cov_dir -1 $fastq1 -r $ref_genome -c $cov_fasta -m $mcclintock_location -p $processors -C $chromosome_names -g $genome
		fi
	else
		if [[ "$remove_intermediates" == "on" ]]
		then
			bash $mcclintock_location/scripts/te_coverage.sh -s $sample -R $referencefolder -o $te_cov_dir -1 $fastq1 -2 $fastq2 -r $ref_genome -c $cov_fasta -m $mcclintock_location -p $processors -i -C $chromosome_names -g $genome
		else
			bash $mcclintock_location/scripts/te_coverage.sh -s $sample -R $referencefolder -o $te_cov_dir -1 $fastq1 -2 $fastq2 -r $ref_genome -c $cov_fasta -m $mcclintock_location -p $processors -C $chromosome_names -g $genome
		fi
	fi
	# unload coverage env and reload conda env for base mcclintock pipeline
	conda deactivate
	conda activate MCCLINTOCK_main
	printf "\nCoverage analysis finished.\n\n" | tee -a /dev/stderr
fi

# If a GFF is supplied then run the analysis using the GFF and TE hierarchy as input
if [[ $inputg ]]
then
	# Copy the te locations file to the run folder
	te_locations_file=${inputg##*/}
	if [[ ! -f $referencefolder/$te_locations_file ]]
	then
		# Copy the gff input via a processing step that creates standard columns and layout for rest of pipeline
		grep -v '^#' $inputg | awk -F'[\t=;]' 'BEGIN {OFS = "\t"}; {printf $1"\t"$2"\t"; for(x=1;x<=NF;x++) if ($x~"ID") printf $(x+1); print "\t"$4,$5,$6,$7,$8,"ID="}' | awk -F'\t' '{print $0$3";Name="$3";Alias="$3}' > $referencefolder/$te_locations_file
	fi
	te_locations=$referencefolder/$te_locations_file

	# Copy the te family file to the run folder
	te_families_file=${inputt##*/}
	if [[ ! -f $referencefolder/$te_families_file ]]
	then
		cp -n $inputt $referencefolder/$te_families_file
	fi
	te_families=$referencefolder/$te_families_file

	# Use the GFF to create input for the rest of the pipeline
	if [[ ! -f $referencefolder/"popoolationte_"$genome".fasta" ]]
	then
		bedtools maskfasta -fi $reference_genome -fo $referencefolder/"popoolationte_"$genome".fasta" -bed $te_locations
	fi
	popoolationte_reference_genome=$referencefolder/"popoolationte_"$genome".fasta"

	# Extract sequence of all reference TE copies if this has not already been done
	# Cut first line if it begins with #
	if [[ ! -f $referencefolder/popool_all_te_seqs.fasta ]]
	then
		# PoPoolationTE always needs the full TE sequences
		bedtools getfasta -name -fi $reference_genome -bed $te_locations -fo $referencefolder/ref_te_seqs.fasta
		cat $consensus_te_seqs $referencefolder/ref_te_seqs.fasta > $referencefolder/popool_te_seqs.fasta
		perl $mcclintock_location/scripts/fixfastalinelength.pl $referencefolder/popool_te_seqs.fasta 80 $referencefolder/popool_all_te_seqs.fasta
		rm $referencefolder/popool_te_seqs.fasta

		if [[ "$addrefcopies" = "on" ]]
		then
			te_seqs=$referencefolder/ref_te_seqs.fasta
		fi
		if [[ "$addconsensus" = "on" ]]
		then
			cat $consensus_te_seqs $te_seqs > $referencefolder/all_te_seqs2.fasta
			te_seqs=$referencefolder/all_te_seqs2.fasta
		fi
		# Use script to fix the line length of reference input to 80 characters (needed for samtools index)
		perl $mcclintock_location/scripts/fixfastalinelength.pl $te_seqs 80 $referencefolder/all_te_seqs.fasta
		rm $reference_genome".fai"
	fi
	all_te_seqs=$referencefolder/all_te_seqs.fasta
	popool_te_seqs=$referencefolder/popool_all_te_seqs.fasta

	# The pipeline functions most comprehensively (i.e. dealing with insertions with no copies in the reference genome) if
	# the sequences of TEs are added to the end of the genome and reflected in the annotation
	if [[ ! -f $referencefolder/full_reference.fasta ]]
	then
		if [[ "$addconsensus" = "on" || "$addrefcopies" = "on" ]]
		then
			cat $reference_genome $all_te_seqs > $referencefolder/full_reference.fasta
			cp $referencefolder/full_reference.fasta $referencefolder/$genome".fasta"
			reference_genome=$referencefolder/$genome".fasta"
		fi
	fi

	# PoPoolationTE always needs the full combination reference
	if [[ ! -f $referencefolder/"popoolationte_full_"$genome".fasta" ]]
	then
			cat $popoolationte_reference_genome $popool_te_seqs > $referencefolder/"popoolationte_full_"$genome".fasta"
	fi
	popoolationte_reference_genome=$referencefolder/"popoolationte_full_"$genome".fasta"

	# Add the locations of the sequences of the consensus TEs to the genome annotation
	if [[ ! -f $referencefolder/TE-names ]]
	then
		awk -F">" '/^>/ {print $2"\t"$2}' $consensus_te_seqs > $referencefolder/tmp
		cat $te_families >> $referencefolder/tmp
		cp $te_families $referencefolder/"popool_"$te_families_file
		cp $te_locations $referencefolder/"popool_"$te_locations_file
		if [[ "$addconsensus" = "on" || "$addrefcopies" = "on" ]]
		then
			awk -F">" '/^>/ {if (seqlen){print seqlen}; printf $2"\t" ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' $all_te_seqs > $referencefolder/TE-lengths
			while read TE length
			do
				echo -e "$TE\treannotate\ttransposable_element\t1\t$length\t.\t+\t.\tID=instance$TE;Name=instance$TE;Alias=instance$TE" >> $te_locations
				awk -vTE=$TE '{ if(TE==$1) print "instance"TE"\t"$2; }' $referencefolder/tmp >> $te_families
			done < $referencefolder/TE-lengths
		fi

		# PoPoolationTE always needs the full family file and annotation
		awk -F">" '/^>/ {print $2}' $consensus_te_seqs > $referencefolder/TE-names
		awk '{ print $1"\t"$1 }' $referencefolder/TE-names >> $referencefolder/"popool_"$te_families_file
		rm $referencefolder/tmp
	fi
	popool_te_locations=$referencefolder/"popool_"$te_locations_file
	popool_te_families=$referencefolder/"popool_"$te_families_file

# The GFF input is optional, if it is not supplied then RepeatMasker is run to generate the necessary inputs
else
	if [[ ! -f $reference_genome".masked" || ! -f $reference_genome".out.gff" ]]
	then
		# Run RepeatMasker on the genome using the TE database to generate gff annotation
		RepeatMasker -pa $processors -lib $consensus_te_seqs -s -gff -nolow -no_is $reference_genome
		# RepeatMasker appears to override the custom database names during the ProcessRepeats step so this changes them back for
		# Drosophila, more rules like this may be needed for other reference genomes
		sed "s/McClintock-int/McClintock/g" $reference_genome".out.gff" > $referencefolder/tmp
		sed "s/POGON1/pogo/g" $referencefolder/tmp > $reference_genome".out.gff"
		perl $mcclintock_location/scripts/fixfastalinelength.pl $reference_genome".masked" 80 $reference_genome".masked2"
		mv $reference_genome".masked2" $reference_genome".masked"
	fi
	popoolationte_reference_genome=$reference_genome".masked"
	te_locations=$reference_genome".out.gff"

	# Run the perl script to create a hierarchy file that corresponds to the RepeatMasker GFF file.
	# (The RepeatMasker file is edited and renamed ..._ID in the process)
	if [[ ! -f $referencefolder/hierarchy.tsv ]]
	then
		perl $mcclintock_location/scripts/hierarchyfromrepeatmasked.pl $te_locations $consensus_te_seqs $referencefolder/hierarchy.tsv
		mv $te_locations"_ID" $te_locations
	fi
	te_families=$referencefolder/hierarchy.tsv
	consensus_te_seqs=$consensus_te_seqs"_ID"

	# Extract sequence of all reference TE copies if this has not already been done
	# Cut first line if it begins with #
	if [[ ! -f $referencefolder/reference_te_seqs.fasta ]]
	then
		grep -v '^#' $te_locations | awk -F'[\t=;]' 'BEGIN {OFS = "\t"}; {printf $1"\t"$2"\t"; for(x=1;x<=NF;x++) if ($x~"ID") printf $(x+1); print "\t"$4,$5,$6,$7,$8,"ID="}' | awk -F'\t' '{print $0$3";Name="$3";Alias="$3}' > edited.gff
		mv edited.gff $te_locations
		bedtools getfasta -name -fi $reference_genome -bed $te_locations -fo $referencefolder/reference_te_seqs.fasta
		# Use script to fix the line length of reference input to 80 characters (needed for samtools index)
		perl $mcclintock_location/scripts/fixfastalinelength.pl $referencefolder/reference_te_seqs.fasta 80 $referencefolder/reference_te_seqs2.fasta
		mv $referencefolder/reference_te_seqs2.fasta $referencefolder/reference_te_seqs.fasta
		rm $reference_genome".fai"
	fi
	reference_te_seqs=$referencefolder/reference_te_seqs.fasta

	if [[ ! -f $referencefolder/popool_te_seqs.fasta ]]
	then
		if [[ "$addconsensus" = "on" ]]
		then
			cat $consensus_te_seqs > $referencefolder/all_te_seqs.fasta
		fi
		if [[ "$addrefcopies" = "on" ]]
		then
			cat $reference_te_seqs >> $referencefolder/all_te_seqs.fasta
		fi
		# PoPoolationTE always needs the full TE sequences
		cat $consensus_te_seqs $reference_te_seqs > $referencefolder/popool_te_seqs.fasta
	fi
	all_te_seqs=$referencefolder/all_te_seqs.fasta
	popool_te_seqs=$referencefolder/popool_te_seqs.fasta

	# Add the locations of the sequences of the consensus TEs to the genome annotation
	if [[ ! -f $referencefolder/TE-lengths ]]
	then
		awk -F">" '/^>/ {print $2"\t"$2}' $consensus_te_seqs > $referencefolder/tmp
		cat $te_families >> $referencefolder/tmp
		cp $te_families $referencefolder/popool_hierarchy.tsv
		cp $te_locations $referencefolder/popool_te_locations.gff
		if [[ "$addconsensus" = "on" ||  "$addrefcopies" = "on" ]]
		then
			awk -F">" '/^>/ {if (seqlen){print seqlen}; printf $2"\t" ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' $all_te_seqs > $referencefolder/TE-lengths
			while read TE length
			do
				echo -e "$TE\treannotate\ttransposable_element\t1\t$length\t.\t+\t.\tID=instance$TE;Name=instance$TE;Alias=instance$TE" >> $te_locations
				awk -vTE=$TE '{ if(TE==$1) print "instance"TE"\t"$2; }' $referencefolder/tmp >> $te_families
			done < $referencefolder/TE-lengths
		fi

		# PoPoolationTE always needs the full family file and annotation
		awk -F">" '/^>/ {print $2}' $consensus_te_seqs > $referencefolder/TE-names
		awk '{ print $1"\t"$1 }' $referencefolder/TE-names >> $referencefolder/popool_hierarchy.tsv
		rm $referencefolder/tmp
	fi
	popool_te_locations=$referencefolder/popool_te_locations.gff
	popool_te_families=$referencefolder/popool_hierarchy.tsv

	# The pipeline functions most comprehensively (i.e. dealing with insertions with no copies in the reference genome) if
	# the sequences of TEs are added to the end of the genome and reflected in the annotation
	if [[ ! -f $referencefolder/full_reference.fasta ]]
	then
		if [[ "$addconsensus" = "on" ||  "$addrefcopies" = "on" ]]
		then
			cat $reference_genome $all_te_seqs > $referencefolder/full_reference.fasta
			cp $referencefolder/full_reference.fasta $referencefolder/$genome".fasta"
			reference_genome=$referencefolder/$genome".fasta"
		fi
	fi

	# PoPoolationTE always needs the full combination reference
	if [[ ! -f $referencefolder/"popoolationte_full_"$genome".fasta" ]]
	then
		cat $popoolationte_reference_genome $popool_te_seqs > $referencefolder/"popoolationte_full_"$genome".fasta"
	fi
	popoolationte_reference_genome=$referencefolder/"popoolationte_full_"$genome".fasta"
fi

# If FastQC is installed then launch FastQC on the input fastqs
location=`which fastqc`
if [[ -z "$location" ]]
then
	printf "\nFastQC not installed, skipping input quality analysis...\n\n" | tee -a /dev/stderr
else
	printf "\nPerforming FastQC analysis...\n\n" | tee -a /dev/stderr
	mkdir $samplefolder/results/summary/fastqc_analysis
	if [[ $single_end == "true" ]]
	then
		fastqc -t $processors $fastq1 -o $samplefolder/results/summary/fastqc_analysis
	else
		fastqc -t $processors $fastq1 $fastq2 -o $samplefolder/results/summary/fastqc_analysis
	fi
fi

# Create indexes of reference genome if not already made for this genome
if [[ ! -f $reference_genome".fai" ]]
then
	samtools faidx $reference_genome
	samtools faidx $popoolationte_reference_genome
fi
if [[ ! -f $reference_genome".bwt" ]]
then
	bwa index $reference_genome
	bwa index $popoolationte_reference_genome
fi

# Create bed file of reference TE locations
if [[ ! -f $referencefolder/reference_TE_locations.bed ]]
then
	awk -F["\t"\;=] '{print $1"\t"$4-1"\t"$5"\t"$10"\t.\t"$7}' $te_locations > $referencefolder/reference_TE_locations.bed
fi
bed_te_locations_file=$referencefolder/reference_TE_locations.bed


# Adjust hierachy levels for TE-locate and TEMP output
telocate_te_locations=${te_locations%.*}
telocate_te_locations=$telocate_te_locations"_HL.gff"
if [[ ! -f $telocate_te_locations ]]
then
	perl $mcclintock_location/TE-locate/TE_hierarchy.pl $te_locations $te_families Alias
fi

# Allow case insensitivity for method names
shopt -s nocasematch
if [[ $methods == *TE-locate* || $methods == *TElocate* || $methods == *RetroSeq* || $methods == *TEMP* ]]
then
	shopt -u nocasematch
	# Create sam files for input
	printf "\nCreating sam alignment...\n\n" | tee -a /dev/stderr

	if [[ "$save_tag" != "on" ]]
	then
		bwa mem -t $processors -v 0 -R '@RG\tID:'$sample'\tSM:'$sample $reference_genome $fastq1 $fastq2 > $samplefolder/sam/$sample.sam
	else
		bwa mem -C -t $processors -v 0 -R '@RG\tID:'$sample'\tSM:'$sample $reference_genome $fastq1 $fastq2 > $samplefolder/sam/$sample.sam
	fi
	sam=$samplefolder/sam/$sample.sam

	# Allow case insensitivity for method names
	shopt -s nocasematch
	if [[ $methods == *TE-locate* || $methods == *TElocate* || $methods == *TEMP* ]]
	then
		shopt -u nocasematch
		# Calculate the median insert size of the sample
		printf "\nCalculating median insert size...\n\n" | tee -a /dev/stderr
		median_insertsize=`cut -f9 $sam | sort -S $memory"G" -n | awk '{if ($1 > 0) ins[reads++]=$1; } END { print ins[int(reads/2)]; }'`
		printf "\nMedian insert size = $median_insertsize\n\n" | tee -a /dev/stderr
		mkdir -p $samplefolder/results/summary
		echo $median_insertsize > $samplefolder/results/summary/median_insertsize
	fi

	# Allow case insensitivity for method names
	shopt -s nocasematch
	if [[ $methods == *RetroSeq* || $methods == *TEMP* ]]
	then
		shopt -u nocasematch
		# Create bam files for input
		printf "\nCreating bam alignment files...\n\n" | tee -a /dev/stderr
		samtools view -Sb -t $reference_genome".fai" $sam | samtools sort - $samplefolder/bam/$sample
		bam=$samplefolder/bam/$sample.bam
		samtools index $bam

		# Get stats of bam file from samtools
		samtools flagstat $bam > $samplefolder/results/summary/bwamem_bamstats.txt
	fi

	# Allow case insensitivity for method names
	shopt -s nocasematch
	if [[ $methods == *TE-locate* || $methods == *TElocate* ]]
	then
		shopt -u nocasematch
		# Sort the sam file lexically for TE-locate
		printf "\nSorting sam alignment...\n\n" | tee -a /dev/stderr
		sort -S $memory"G" --temporary-directory=$samplefolder/sam/ $samplefolder/sam/$sample.sam > $samplefolder/sam/sorted$sample.sam
		rm $samplefolder/sam/$sample.sam
		mv $samplefolder/sam/sorted$sample.sam $samplefolder/sam/$sample.sam
		sam=$samplefolder/sam/$sample.sam
		sam_folder=$samplefolder/sam
	else
		# If TE-locate isn't going to be launched then the SAM file is no longer needed
		if [[ "$remove_intermediates" = "on" ]]
		then
			rm -r $sam_folder
		fi
	fi
fi

# Allow case insensitivity for method names
shopt -s nocasematch
if [[ $methods == *TE-locate* || $methods == *TElocate* ]]
then
	shopt -u nocasematch
	################################## Run TE-locate ##################################

	printf "\nRunning TE-locate pipeline...\n\n" | tee -a /dev/stderr

	# If there are not at least five chromosomes in the reference genome then TE-locate will fail
	telocate_reference_genome=$referencefolder/$genome"TElocate.fasta"
	if [[ ! -f $telocate_reference_genome ]]
	then
		cp $reference_genome $telocate_reference_genome
		# If there are less than 5 chromosomes, top them up to 5 with short false ones
		chromosomes=`wc -l $chromosome_names | awk '{print $1}'`
		if [[ $chromosomes < 5 ]]
		then
			needed=$((5 - $chromosomes))
			for i in `seq 1 $needed`
			do
				echo ">fixforTElocate"$i >> $telocate_reference_genome
				echo "ACGT" >> $telocate_reference_genome
			done
		fi
	fi

	cd $mcclintock_location/TE-locate

	# set up conda env for telocate
	conda deactivate
	conda activate MCCLINTOCK_telocate
	bash runtelocate.sh $sam_folder $telocate_reference_genome $telocate_te_locations $memory $sample $median_insertsize $samplefolder
	# unload telocate conda env and reload env for base mcclintock
	conda deactivate
	conda activate MCCLINTOCK_main

	# Save the original result file and the bed files filtered by mcclintock
	mv $samplefolder/TE-locate/$sample"_telocate"* $samplefolder/results/
	mkdir $samplefolder/results/originalmethodresults/TE-locate
	cp $samplefolder/TE-locate/*.info $samplefolder/results/originalmethodresults/TE-locate

	# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
	# i.e. leave any reusable species data behind.
	if [[ "$remove_intermediates" = "on" ]]
	then
		printf "\nRemoving TE-locate intermediate files\n\n" | tee -a /dev/stderr
		rm -r $sam_folder
		rm -r $samplefolder/TE-locate
	fi
	cd $mcclintock_location/
fi

# Allow case insensitivity for method names
shopt -s nocasematch
if [[ $methods == *RetroSeq* ]]
then
	shopt -u nocasematch
	################################## Run RetroSeq ##################################

	printf "\nRunning RetroSeq pipeline...\n\n" | tee -a /dev/stderr
	
	cd $mcclintock_location/RetroSeq

	# set up conda env for retroseq
	conda deactivate
	conda activate MCCLINTOCK_retroseq
	bash runretroseq.sh $consensus_te_seqs $bam $reference_genome $bed_te_locations_file $te_families $samplefolder
	# unload retroseq conda env and reload env for base mcclintock
	conda deactivate
	conda activate MCCLINTOCK_main

	# Save the original result file and the bed files filtered by mcclintock
	mv $samplefolder/RetroSeq/$sample"_retroseq"* $samplefolder/results/
	mkdir $samplefolder/results/originalmethodresults/RetroSeq
	cp $samplefolder/RetroSeq/$sample".calling.PE.vcf" $samplefolder/results/originalmethodresults/RetroSeq

	# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
	# i.e. leave any reusable species data behind.
	if [[ "$remove_intermediates" = "on" ]]
	then
		# If the save bam option is specified then override the command to delete the bam files.
		if [[ "$save_bam" != "on" ]]
		then
			shopt -s nocasematch
			if [[ $methods != *TEMP* ]]
			then
				shopt -u nocasematch
				rm -r $samplefolder/bam
			fi
		fi
		printf "\nRemoving RetroSeq intermediate files\n\n" | tee -a /dev/stderr
		rm -r $samplefolder/RetroSeq
	fi
	cd $mcclintock_location/
fi

# Allow case insensitivity for method names
shopt -s nocasematch
if [[ $methods == *TEMP* ]]
then
	shopt -u nocasematch
	################################## Run TEMP ##################################

	printf "\nRunning TEMP pipeline...\n\n" | tee -a /dev/stderr

	if [[ ! -f $referencefolder/$genome.2bit ]]
	then
		faToTwoBit $reference_genome $referencefolder/$genome.2bit
	fi
	twobitreference=$referencefolder/$genome.2bit

	cd $mcclintock_location/TEMP

	# set up conda env for temp
	conda deactivate
	conda activate MCCLINTOCK_temp
	bash runtemp.sh $bam $twobitreference $consensus_te_seqs $bed_te_locations_file $te_families $median_insertsize $sample $processors $samplefolder $telocate_te_locations
	# unload temp conda env and reload env for base mcclintock
	conda deactivate
	conda activate MCCLINTOCK_main

	# Save the original result file and the bed files filtered by mcclintock
	mv $samplefolder/TEMP/$sample"_temp"* $samplefolder/results/
	mkdir $samplefolder/results/originalmethodresults/TEMP
	cp $samplefolder/TEMP/$sample".insertion.refined.bp.summary" $samplefolder/results/originalmethodresults/TEMP
	cp $samplefolder/TEMP/$sample".absence.refined.bp.summary" $samplefolder/results/originalmethodresults/TEMP

	# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
	# i.e. leave any reusable species data behind.
	if [[ "$remove_intermediates" = "on" ]]
	then
		# If the save bam option is specified then override the command to delete the bam files.
		if [[ "$save_bam" != "on" ]]
		then
			rm -r $samplefolder/bam
		fi
		printf "\nRemoving TEMP intermediate files\n\n" | tee -a /dev/stderr
		rm -r $samplefolder/TEMP
	fi
	cd $mcclintock_location/
fi

# Allow case insensitivity for method names
shopt -s nocasematch
if [[ $methods == *RelocaTE* ]]
then
	shopt -u nocasematch
	################################## Run RelocaTE ##################################

	printf "\nRunning RelocaTE pipeline...\n\n" | tee -a /dev/stderr

	# Add TSD lengths to consensus TE sequences
	if [[ ! -f $referencefolder/relocate_te_seqs.fasta ]]
	then
		awk '{if (/>/) print $0" TSD=UNK"; else print $0}' $consensus_te_seqs > $referencefolder/relocate_te_seqs.fasta
	fi
	relocate_te_seqs=$referencefolder/relocate_te_seqs.fasta

	# Create general gff file to allow reference TE detection in RelocaTE
	if [[ ! -f $referencefolder/relocate_te_locations.gff ]]
	then
		awk 'FNR==NR{array[$1]=$2;next}{print $1,$2,array[$3],$4,$5,$6,$7,$8,$9}' FS='\t' OFS='\t' $te_families $te_locations > $referencefolder/relocate_te_locations.gff
	fi
	relocate_te_locations=$referencefolder/relocate_te_locations.gff

	# set up conda env for relocate
	cd $mcclintock_location/RelocaTE
	conda deactivate
	conda activate MCCLINTOCK_relocate
	bash runrelocate.sh $relocate_te_seqs $reference_genome $samplefolder/reads $sample $relocate_te_locations $samplefolder $single_end
	# unload relocate conda env and reload env for base mcclintock
	conda deactivate
	conda activate MCCLINTOCK_main

	# Save the original result file and the bed files filtered by mcclintock
	mv $samplefolder/RelocaTE/$sample"_relocate"* $samplefolder/results/
	mkdir $samplefolder/results/originalmethodresults/RelocaTE
	cp -r $samplefolder/RelocaTE/*/results $samplefolder/results/originalmethodresults/RelocaTE

	# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
	# i.e. leave any reusable species data behind.
	if [[ "$remove_intermediates" = "on" ]]
	then
		printf "\nRemoving RelocaTE intermediate files\n\n" | tee -a /dev/stderr
		rm -r $samplefolder/RelocaTE
	fi
	cd $mcclintock_location/
fi

# Allow case insensitivity for method names
shopt -s nocasematch
if [[ $methods == *ngs_te_mapper* || $methods == *ngstemapper* ]]
then
	shopt -u nocasematch
	################################## Run ngs_te_mapper ##################################

	printf "\nRunning ngs_te_mapper pipeline...\n\n" | tee -a /dev/stderr

	cd $mcclintock_location/ngs_te_mapper

	# If the data is only single ended then pass the message "false" to ngs_te_mapper
	if [[ $single_end == "true" ]]
	then
		fastq2="false"
	fi

	# set up conda env for ngs_te_mapper
	conda deactivate
	conda activate MCCLINTOCK_ngstemapper
	bash runngstemapper.sh $consensus_te_seqs $reference_genome $sample $fastq1 $fastq2 $samplefolder $processors
	# unload ngs_te_mapper conda env and reload env for base mcclintock
	conda deactivate
	conda activate MCCLINTOCK_main

	# Save the original result file and the bed file filtered by mcclintock
	mv $samplefolder/ngs_te_mapper/$sample"_ngs_te_mapper_nonredundant.bed" $samplefolder/results/
	mkdir $samplefolder/results/originalmethodresults/ngs_te_mapper
	cp $samplefolder/ngs_te_mapper/bed_tsd/*.bed $samplefolder/results/originalmethodresults/ngs_te_mapper

	# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
	# i.e. leave any reusable species data behind.
	if [[ "$remove_intermediates" = "on" ]]
	then
		printf "\nRemoving ngs_te_mapper intermediate files\n\n" | tee -a /dev/stderr
		rm -r $samplefolder/ngs_te_mapper
	fi
	cd $mcclintock_location/
fi

# Allow case insensitivity for method names
shopt -s nocasematch
if [[ $methods == *popoolationte* ]]
then
	shopt -u nocasematch
	################################## Run PoPoolationTE ##################################

	printf "\nRunning PoPoolationTE pipeline...\n\n" | tee -a /dev/stderr

	# Create te_hierachy
	if [[ ! -f $referencefolder/te_hierarchy ]]
	then
		printf "insert\tid\tfamily\tsuperfamily\tsuborder\torder\tclass\tproblem\n" > $referencefolder/te_hierarchy
		awk '{printf $0"\t"$2"\t"$2"\tna\tna\tna\t0\n"}' $popool_te_families >> $referencefolder/te_hierarchy
	fi
	te_hierarchy=$referencefolder/te_hierarchy

	cd $mcclintock_location/PoPoolationTE
	
	# set up conda env for populationte
	conda deactivate
	conda activate MCCLINTOCK_popoolationte
	bash runpopoolationte.sh $popoolationte_reference_genome $te_hierarchy $fastq1 $fastq2 $popool_te_locations $processors $samplefolder
	# unload populationte conda env and reload env for base mcclintock
	conda deactivate
	conda activate MCCLINTOCK_main

	# Save the original result file and the bed files filtered by mcclintock
	mv $samplefolder/PoPoolationTE/$sample"_popoolationte"* $samplefolder/results/
	mkdir $samplefolder/results/originalmethodresults/PoPoolationTE
	cp $samplefolder/PoPoolationTE/te-poly-filtered.txt $samplefolder/results/originalmethodresults/PoPoolationTE

	# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
	# i.e. leave any reusable species data behind.
	if [[ "$remove_intermediates" = "on" ]]
	then
		printf "\nRemoving PoPoolationTE intermediate files\n\n" | tee -a /dev/stderr
		rm -r $samplefolder/PoPoolationTE
	fi
	cd $mcclintock_location/
fi

# Generate final report
## McClintock version and invoke command
report=$samplefolder/results/summary/summary_report.txt
echo "McClintock version: `git rev-parse HEAD`" > $report
echo "McClintock was called as follows:" >> $report
invocation="$(printf %q "$BASH_SOURCE")$((($#)) && printf ' %q' "$@")"
echo "$invocation" >> $report

## read length
if [[ -n "$location" ]]
then
	fastqc_dir=$samplefolder/results/summary/fastqc_analysis
	cd $fastqc_dir
	unzip '*.zip'
	if [[ $single_end != "true" ]]
	then
		r1_fastqc=`ls -d -1 $fastqc_dir/*_1*/fastqc_data.txt`
		r2_fastqc=`ls -d -1 $fastqc_dir/*_2*/fastqc_data.txt`
		r1_len=`grep "Sequence length" $r1_fastqc | sed "s/Sequence length/read1 sequence length:/g"`
		r2_len=`grep "Sequence length" $r2_fastqc | sed "s/Sequence length/read2 sequence length:/g"`
		echo "$r1_len" >> $report
		echo "$r2_len" >> $report
		rm -rf "$fastqc_dir/$sample"_1_fastqc
		rm -rf "$fastqc_dir/$sample"_2_fastqc
	else
		r1_fastqc=`ls -d -1 $fastqc_dir/*_1*/fastqc_data.txt`
		r1_len=`grep "Sequence length" $r1_fastqc | sed "s/Sequence length/sequence length:/g"`
		echo "$r1_len" >> $report
		rm -rf "$fastqc_dir/$sample"_fastqc
	fi
fi

## read count
bamstat=$samplefolder/results/summary/bwamem_bamstats.txt
if [[ -e "$bamstat" ]]
then
	r1_count=`grep "read1" $bamstat | sed 's/\s.*//g'`
	r2_count=`grep "read2" $bamstat | sed 's/\s.*//g'`
	if [[ -z "$r2_count" ]]
	then
		echo "Number of reads: $r1_count" >> $report
	else
		echo "Number of reads in read1: $r1_count" >> $report
		echo "Number of reads in read2: $r2_count" >> $report
	fi
fi

## median insert size
median_insertsize=$samplefolder/results/summary/median_insertsize
if [[ -e "$median_insertsize" ]]
then
	median_ins=`cat $median_insertsize`
	echo "median insert size: $median_ins" >> $report
fi

## average genome coverage
if [[ -e "$bam" ]]
then
	if [[ ! -f $referencefolder/"$genome.fasta.out.complement.all.bed" ]]
	then
		cat "$reference_genome".fai | cut -f1,2 | sort -k 1,1 -k2,2n > $referencefolder/"$genome.sorted.all.genome"
		bedtools slop -i $bed_te_locations_file -g $referencefolder/"$genome.sorted.all.genome" -l 1 -r 0 | sort -k 1,1 -k2,2n > $referencefolder/"$genome.fasta.out.zero.all.bed"
		bedtools complement -i $referencefolder/"$genome.fasta.out.zero.all.bed" -g $referencefolder/"$genome.sorted.all.genome" | grep -v "\b0\s*0\b" | grep -v "\s-[0-9]*\b" > $referencefolder/"$genome.fasta.out.complement.all.bed"
	fi
	bed_nonte=$referencefolder/"$genome.fasta.out.complement.all.bed"
	genome_avg_depth=`samtools depth -b  $bed_nonte $bam | awk '{ total += $3 } END { print total/NR }'`
	echo "average genome coverage: $genome_avg_depth" >> $report
fi

# per sample number of insertions (ref and non-ref)
ngs_out="$samplefolder/results/$sample"_ngs_te_mapper_nonredundant.bed
if [[ -e "$ngs_out" ]]
then
	te_all=`tail -n +2 $ngs_out | wc -l`
	te_non_ref=`grep -c "non-ref" $ngs_out`
	te_ref=`tail -n +2 $ngs_out | grep -v -c "non-ref"`
	echo "TEs detected by ngs_te_mapper (all): $te_all" >> $report
	echo "TEs detected by ngs_te_mapper (non-reference): $te_non_ref" >> $report
	echo "TEs detected by ngs_te_mapper (reference): $te_ref" >> $report
fi

relocate_out="$samplefolder/results/$sample"_relocate_nonredundant.bed
if [[ -e "$relocate_out" ]]
then
	te_all=`tail -n +2 $relocate_out | wc -l`
	te_non_ref=`grep -c "non-ref" $relocate_out`
	te_ref=`tail -n +2 $relocate_out | grep -v -c "non-ref"`
	echo "TEs detected by RelocaTE (all): $te_all" >> $report
	echo "TEs detected by RelocaTE (non-reference): $te_non_ref" >> $report
	echo "TEs detected by RelocaTE (reference): $te_ref" >> $report
fi

temp_out="$samplefolder/results/$sample"_temp_nonredundant.bed
if [[ -e "$temp_out" ]]
then
	te_all=`tail -n +2 $temp_out | wc -l`
	te_non_ref=`grep -c "non-ref" $temp_out`
	te_ref=`tail -n +2 $temp_out | grep -v -c "non-ref"`
	echo "TEs detected by TEMP (all): $te_all" >> $report
	echo "TEs detected by TEMP (non-reference): $te_non_ref" >> $report
	echo "TEs detected by TEMP (reference): $te_ref" >> $report
fi

popoolationte_out="$samplefolder/results/$sample"_popoolationte_nonredundant.bed
if [[ -e "$popoolationte_out" ]]
then
	te_all=`tail -n +2 $popoolationte_out | wc -l`
	te_non_ref=`grep -c "non-ref" $popoolationte_out`
	te_ref=`tail -n +2 $popoolationte_out | grep -v -c "non-ref"`
	echo "TEs detected by PoPoolationTE (all): $te_all" >> $report
	echo "TEs detected by PoPoolationTE (non-reference): $te_non_ref" >> $report
	echo "TEs detected by PoPoolationTE (reference): $te_ref" >> $report
fi

retroseq_out="$samplefolder/results/$sample"_retroseq_nonredundant.bed
if [[ -e "$retroseq_out" ]]
then
	te_all=`tail -n +2 $retroseq_out | wc -l`
	te_non_ref=`grep -c "non-ref" $retroseq_out`
	te_ref=`tail -n +2 $retroseq_out | grep -v -c "non-ref"`
	echo "TEs detected by RetroSeq (all): $te_all" >> $report
	echo "TEs detected by RetroSeq (non-reference): $te_non_ref" >> $report
	echo "TEs detected by RetroSeq (reference): $te_ref" >> $report
fi

telocate_out="$samplefolder/results/$sample"_telocate_nonredundant.bed
if [[ -e "$telocate_out" ]]
then
	te_all=`tail -n +2 $telocate_out | wc -l`
	te_non_ref=`grep -c "non-ref" $telocate_out`
	te_ref=`tail -n +2 $telocate_out | grep -v -c "non-ref"`
	echo "TEs detected by TE-locate (all): $te_all" >> $report
	echo "TEs detected by TE-locate (non-reference): $te_non_ref" >> $report
	echo "TEs detected by TE-locate (reference): $te_ref" >> $report
fi

# create summary table for family and method based TE detection
if [[ $methods != "" ]]
then
	bed_dir="$samplefolder/results"
	te_name=$referencefolder/TE-names
	te_summary_out="$samplefolder/results/summary/te_summary_table.csv"
	perl $mcclintock_location/scripts/te_summary_table.pl -d $bed_dir -t $te_name > $te_summary_out
fi
#########################################################################################

# If a user has used an altered genome then insertions in false chromosomes may exist
# These can be because of a real nested insertion or self similarity
# These results are removed to prevent errors in later analysis (e.g. genome browsing) but are saved in a folder
if [[ $methods != "" ]]
then
	if [[ "$addconsensus" = "on" || "$addrefcopies" = "on" ]]
	then
		mkdir $samplefolder/results/non-ref_chromosome_results
		for result_file in $samplefolder/results/*.bed
		do
			result_file_name=`basename $result_file`
			grep  -w -v -F -f $chromosome_names $result_file > $samplefolder/results/non-ref_chromosome_results/$result_file_name
			head -1 $result_file > $samplefolder/results/reference_chr_results
			grep  -w -F -f $chromosome_names $result_file >> $samplefolder/results/reference_chr_results
			mv $samplefolder/results/reference_chr_results $result_file
		done
	fi
fi
# If cleanup intermediate files is specified then delete all intermediate files specific to the sample
# i.e. leave any reusable species data behind.
if [[ "$remove_intermediates" = "on" ]]
then
	printf "\nRemoving McClintock intermediate files\n\n" | tee -a /dev/stderr
	rm -r $samplefolder/reads
fi

# unload conda env
conda deactivate

printf "\nPipeline Complete\n\n" | tee -a /dev/stderr
