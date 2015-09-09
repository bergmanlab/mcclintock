#!/bin/bash

# Download all of the software methods.
echo "Downloading pipelines..."
wget --no-check-certificate https://popoolationte.googlecode.com/files/popoolationte_1.02.zip -O PoPoolationTE.zip 
wget --no-check-certificate https://github.com/bergmanlab/ngs_te_mapper/archive/79ef861f1d52cdd08eb2d51f145223fad0b2363c.zip -O ngs_te_mapper.zip
wget --no-check-certificate https://github.com/tk2/RetroSeq/archive/700d4f76a3b996686652866f2b81fefc6f0241e0.zip -O RetroSeq.zip
wget --no-check-certificate https://github.com/srobb1/RelocaTE/archive/ce3a2066e15f5c14e2887fdf8dce0485e1750e5b.zip -O RelocaTE.zip
wget http://sourceforge.net/projects/te-locate/files/TE-locate.tar/download -O TE-locate.tar
wget --no-check-certificate https://github.com/JialiUMassWengLab/TEMP/archive/d2500b904e2020d6a1075347b398525ede5feae1.zip -O TEMP.zip

# Extract software and format folder layout.
echo "Extracting pipelines..." 
unzip PoPoolationTE.zip
mv popoolationte PoPoolationTE
rm PoPoolationTE.zip
unzip ngs_te_mapper.zip
mv ngs_te_mapper-79ef861f1d52cdd08eb2d51f145223fad0b2363c ngs_te_mapper
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

# Apply edits to software and custom run scripts.
echo "Copying run scripts..."
cp scripts/runpopoolationte.sh PoPoolationTE
patch PoPoolationTE/Modules/TEInsertUtility.pm < scripts/TEInsertUtility.patch
patch PoPoolationTE/Modules/TEInsert.pm < scripts/TEInsert.patch
patch PoPoolationTE/samro.pl < scripts/samro.patch

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
patch TEMP/scripts/TEMP_Absense.sh < TEMP_Absence.patch

cp scripts/mcclintock.sh .

# Change #! to make software more compatible with different environments
sed -i 's/#!\/usr\/bin\/perl -w/#!\/usr\/bin\/env perl/' RelocaTE/scripts/*.pl
# Remove deprecated use of defined from RelocaTE to reduce errors
for f in RelocaTE/scripts/*.pl
do
	sed -i 's/defined @/@/g' $f
done

# Check dependencies 
echo "Testing dependencies..."
dependencies="R RepeatMasker bedtools samtools bcftools bwa exonerate bowtie blat twoBitToFa java perl"

for dependency in $dependencies
do 
	location=`which $dependency`
	if test -z "$location"
	then
		echo "$dependency... NOT FOUND"
	else 
		echo "$dependency... FOUND"
		if [[ $dependency == "bwa" ]]
		then
			version=`bwa 2>&1 >/dev/null | grep "Version:"`
			if [[ $version != "Version: 0.7.4-r385" ]]
			then
				echo "CAUTION McClintock requires version 0.7.4-r385 to work correctly with all methods"
				echo "You currently have $version installed."
			fi
		fi
	fi
done

# Test BioPerl
bioperltest=$(perl -e 'use Bio::Seq' 2>&1)
if [[ -n "$bioperltest" ]]
then
	echo "BioPerl... NOT FOUND"
else  
	echo "BioPerl... FOUND"
fi

# Test for fastqc - not required but recommended
echo "Testing optional dependencies..."
location=`which fastqc`
if test -z "$location"
then
	echo "fastqc... NOT FOUND (Not required but recommended)"
else
	echo "fastqc... FOUND"
fi
