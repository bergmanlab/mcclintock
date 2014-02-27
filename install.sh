#!/bin/bash

# Download all of the software methods.
echo "Downloading pipelines..."
wget --no-check-certificate https://popoolationte.googlecode.com/files/popoolationte_1.02.zip -O PoPoolationTE.zip 
wget --no-check-certificate https://github.com/bergmanlab/ngs_te_mapper/archive/52cfa4d0eefd9aca886f214629b33cfbcb6ebf45.zip -O ngs_te_mapper.zip
wget --no-check-certificate https://github.com/tk2/RetroSeq/archive/700d4f76a3b996686652866f2b81fefc6f0241e0.zip -O RetroSeq.zip
wget --no-check-certificate https://github.com/srobb1/RelocaTE/archive/9b3c89d03c1a8f163d6d2df937f9dc3f02e7e344.zip -O RelocaTE.zip
wget -O TE-locate.tar http://sourceforge.net/projects/te-locate/files/latest/download?source=files

# Extract software and format folder layout.
echo "Extracting pipelines..." 
unzip PoPoolationTE.zip 
rm PoPoolationTE.zip
unzip ngs_te_mapper.zip
mv ngs_te_mapper-52cfa4d0eefd9aca886f214629b33cfbcb6ebf45 ngs_te_mapper
rm ngs_te_mapper.zip
unzip RetroSeq.zip
mv RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0 RetroSeq
rm RetroSeq.zip
unzip RelocaTE.zip
mv RelocaTE-9b3c89d03c1a8f163d6d2df937f9dc3f02e7e344 RelocaTE
rm RelocaTE.zip
mkdir TE-locate
tar -xvf TE-locate.tar -C TE-locate
rm TE-locate.tar

# Apply edits to software and custom run scripts.
echo "Copying run scripts..."
cp scripts/runpopoolationte.sh popoolationte
cp scripts/samro.pl popoolationte
cp scripts/runngstemapper.sh ngs_te_mapper
cp scripts/runretroseq.sh RetroSeq
cp scripts/splitforRetroSeq.pl RetroSeq
cp scripts/runrelocate.sh RelocaTE
cp scripts/runtelocate.sh TE-locate

# Change #! to make software more compatible with different environments
sed -i 's/#!\/usr\/bin\/perl -w/#!\/usr\/bin\/env perl/' RelocaTE/scripts/*.pl

# Check dependencies 
echo "Testing dependencies..."
dependencies="R RepeatMasker bedtools samtools bcftools bwa exonerate bowtie blat perl"

for dependency in $dependencies
do 
	location=`which $dependency`
	if test -z "$location"
	then
		echo "$dependency... NOT FOUND"
	else 
		echo "$dependency... FOUND"
        fi
done

# Test BioPerl
bioperltest=$(perl -e 'use Bio::Seq' 2>&1)
if [ -n "$bioperltest" ]
then
  echo "BioPerl... NOT FOUND"
else  
  echo "BioPerl... FOUND"
fi 
