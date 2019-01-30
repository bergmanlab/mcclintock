#!/bin/bash

# Download all of the software methods.
echo "Downloading pipelines..."
wget --no-check-certificate http://downloads.sourceforge.net/project/popoolationte/popoolationte_1.02.zip -O PoPoolationTE.zip
wget --no-check-certificate https://github.com/bergmanlab/ngs_te_mapper/archive/fb23590200666fe66f1c417c5d5934385cb77ab9.zip -O ngs_te_mapper.zip
wget --no-check-certificate https://github.com/tk2/RetroSeq/archive/700d4f76a3b996686652866f2b81fefc6f0241e0.zip -O RetroSeq.zip
wget --no-check-certificate https://github.com/srobb1/RelocaTE/archive/ce3a2066e15f5c14e2887fdf8dce0485e1750e5b.zip -O RelocaTE.zip
wget --no-check-certificate https://downloads.sourceforge.net/project/te-locate/TE-locate.tar -O TE-locate.tar
wget --no-check-certificate https://github.com/JialiUMassWengLab/TEMP/archive/d2500b904e2020d6a1075347b398525ede5feae1.zip -O TEMP.zip
wget --no-check-certificate https://github.com/bergmanlab/samplot/archive/1de65afd22e88c5cb5122ae638e8ba4cf6f75144.zip -O samplot.zip

# Extract software and format folder layout.
echo "Extracting pipelines..."
unzip PoPoolationTE.zip
mv popoolationte PoPoolationTE
rm PoPoolationTE.zip
unzip ngs_te_mapper.zip
mv ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9 ngs_te_mapper
rm ngs_te_mapper.zip
unzip RetroSeq.zip
mv RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0 RetroSeq
rm RetroSeq.zip
unzip RelocaTE.zip
mv RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b RelocaTE
rm RelocaTE.zip
mkdir TE-locate
tar -xvf TE-locate.tar -C TE-locate
rm TE-locate.tar
unzip TEMP.zip
mv TEMP-d2500b904e2020d6a1075347b398525ede5feae1 TEMP
rm TEMP.zip
unzip samplot.zip
mv samplot-1de65afd22e88c5cb5122ae638e8ba4cf6f75144 samplot
rm samplot.zip

# Apply edits to software and custom run scripts.
echo "Copying run scripts..."
cp scripts/runpopoolationte.sh PoPoolationTE
patch PoPoolationTE/Modules/TEInsertUtility.pm < scripts/TEInsertUtility.patch
patch PoPoolationTE/Modules/TEInsert.pm < scripts/TEInsert.patch
patch PoPoolationTE/samro.pl < scripts/samro.patch
patch PoPoolationTE/identify-te-insertsites.pl < scripts/identify-te-insertsites.patch

cp scripts/runngstemapper.sh ngs_te_mapper

cp scripts/runretroseq.sh RetroSeq
cp scripts/splitforRetroSeq.pl RetroSeq

cp scripts/runrelocate.sh RelocaTE
patch RelocaTE/scripts/relocaTE_insertionFinder.pl < scripts/relocaTE_insertionFinder.patch

cp scripts/runtelocate.sh TE-locate
# These fixes are now built in to TE-locate
#mv scripts/genomes1001_data_processing TE-locate
#cd TE-locate
#jar uf TE_locate_core.jar genomes1001_data_processing/Functions.class
#cd ..

cp scripts/runtemp.sh TEMP
patch TEMP/scripts/TEMP_Absence.sh < scripts/TEMP_Absence.patch

cp scripts/mcclintock.sh .

# Replace PATH for perl in shebang with the one in current working enrivonment
source activate MCCLINTOCK_relocate
perlPath=$(which perl | sed 's/\//\\\//g')
sed -i 's/\/usr\/bin\/perl -w/'"$perlPath"'/' RelocaTE/scripts/*.pl
source deactivate

# Remove deprecated use of defined from RelocaTE to reduce errors
for f in RelocaTE/scripts/*.pl
do
	sed -i 's/defined @/@/g' $f
done

# Check dependencies
check_bioperl() {
	bioperltest=$(perl -e 'use Bio::Seq' 2>&1)
	if [[ -n "$bioperltest" ]]
	then
		echo "BioPerl... NOT FOUND"
	else
		echo "BioPerl... FOUND"
	fi
}

check_bwa_version() {
	version=`bwa 2>&1 >/dev/null | grep "Version:"`
	if [[ $version != "Version: 0.7.4-r385" ]]
	then
		echo "CAUTION McClintock requires version 0.7.4-r385 to work correctly with all methods"
		echo "You currently have $version installed."
	fi
}

check_fastqc() {
	location=`which fastqc`
	if test -z "$location"
	then
		echo "fastqc... NOT FOUND (Not required but recommended)"
	else
		echo "fastqc... FOUND"
	fi
}

check_dependency() {
	dependencies=("$@")
	for dependency in "${dependencies[@]}"
	do
		if [[ $dependency == "bioperl" ]]
		then
			check_bioperl
			continue
		fi
		if [[ $dependency == "fastqc" ]]
		then
			check_fastqc
			continue
		fi
		location=`which $dependency`
		if test -z "$location"
		then
			echo "$dependency... NOT FOUND"
		else
			echo "$dependency... FOUND"
			if [[ $dependency == "bwa" ]]
			then
				check_bwa_version
			fi
		fi
	done
}

echo "Testing dependencies for McClintock main pipeline..."
dependencies_main=("RepeatMasker" "perl" "bwa" "bedtools" "samtools" "fastqc" "faToTwoBit")
source activate MCCLINTOCK_main
check_dependency "${dependencies_main[@]}"
source deactivate

echo "Testing dependencies for McClintock coverage module..."
dependencies_cov=("RepeatMasker" "perl" "bwa" "bedtools" "samtools" "python")
source activate MCCLINTOCK_cov
check_dependency "${dependencies_cov[@]}"
source deactivate

echo "Testing dependencies for McClintock ngs_te_mapper module..."
dependencies_ngstemapper=("bwa" "bedtools" "R")
source activate MCCLINTOCK_ngstemapper
check_dependency "${dependencies_ngstemapper[@]}"
source deactivate

echo "Testing dependencies for McClintock TEMP module..."
dependencies_temp=("perl" "bwa" "bedtools" "samtools" "faToTwoBit" "twoBitToFa")
source activate MCCLINTOCK_temp
check_dependency "${dependencies_temp[@]}"
source deactivate

echo "Testing dependencies for McClintock RelocaTE module..."
dependencies_relocate=("perl" "bedtools" "bowtie" "blat" "samtools")
source activate MCCLINTOCK_relocate
check_dependency "${dependencies_relocate[@]}"
source deactivate

echo "Testing dependencies for McClintock TE-locate module..."
dependencies_telocate=("perl" "bedtools" "bwa" "java")
source activate MCCLINTOCK_telocate
check_dependency "${dependencies_telocate[@]}"
source deactivate

echo "Testing dependencies for McClintock PoPoolationTE module..."
dependencies_popoolationte=("perl" "bedtools" "bwa" "samtools")
source activate MCCLINTOCK_popoolationte
check_dependency "${dependencies_popoolationte[@]}"
source deactivate

echo "Testing dependencies for McClintock RetroSeq module..."
dependencies_retroseq=("perl" "bedtools" "bwa" "samtools" "bcftools" "exonerate")
source activate MCCLINTOCK_retroseq
check_dependency "${dependencies_retroseq[@]}"
source deactivate