Introduction
------
Many methods have been developed to detect transposable element (TE) insertions from whole genome shotgun next-generation sequencing (NGS) data, each of which has different dependencies, run interfaces, and output formats. Here, we have developed a meta-pipeline to download, install and run six available methods for detecting TE insertions in NGS data, which generates output in the UCSC Browser extensible data (BED) format.

Software Methods
------
 * [ngs_te_mapper](https://github.com/bergmanlab/ngs_te_mapper "Click to go to download location") - [Linheiro and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371/journal.pone.0030008 "Click to go to paper location")
 * [PoPoolationTE](https://code.google.com/p/popoolationte/ "Click to go to download location") - [Kofler *et al.* (2012)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002487;jsessionid=2CFC9BF7DEF785D90070915204B5F846 "Click to go to paper location")
 * [RelocaTE](https://github.com/srobb1/RelocaTE "Click to go to download location") - [Robb *et al.* (2013)](http://www.g3journal.org/content/3/6/949.long "Click to go to paper location")
 * [TE-locate](http://zendto.gmi.oeaw.ac.at/pickup.php?claimID= Y3tZVfN5xipYyBDN&claimPasscode=NArXMbTjmkorWjSM&emailAddr=te_locate%40gmx.at "Click to go to download location") - [Platzer *et al.* (2012)](http://www.mdpi.com/2079-7737/1/2/395 "Click to go to paper location")
 * [RetroSeq](https://github.com/tk2/RetroSeq "Click to go to download location") - [Keane *et al.* (2012)](http://bioinformatics.oxfordjournals.org/content/29/3/389.long "Click to go to paper location")
 * [TEMP](https://github.com/JialiUMassWengLab/TEMP "Click to go to download location") - [Zhuang *et al.* (2014)](http://nar.oxfordjournals.org/content/42/11/6826.full "Click to go to paper location")

Software Dependencies
------
All of the software systems must be run on a unix based system with the software dependencies listed per method below. The versions used to run this pipeline are indicated in parentheses and no guarantee is made that it will function using alternate versions.

 * ngs_te_mapper
  * [R](http://cran.r-project.org/) (v.3.0.2)
  * [BWA](http://sourceforge.net/projects/bio-bwa/files/) (v.0.7.2-r351)

 * PoPoolationTE
  * [Perl](http://www.perl.org/get.html) (v.5.14.2)
  * [RepeatMasker](http://www.repeatmasker.org/RMDownload.html) (v.4.0.2)
  * [SAMtools](http://sourceforge.net/projects/samtools/files/) (v.0.1.19-44428cd)
  * [BWA](http://sourceforge.net/projects/bio-bwa/files/) (v.0.7.2-r351)

 * RelocaTE
  * [Perl](http://www.perl.org/get.html) (v.5.14.2)
  * [BioPerl](http://www.bioperl.org/wiki/Getting_BioPerl) (v.1.006901)
  * [SAMtools](http://sourceforge.net/projects/samtools/files/) (v.0.1.19-44428cd)
  * [Blat](http://users.soe.ucsc.edu/~kent/src/) (v.35x1)
  * [Bowtie](http://bowtie-bio.sourceforge.net/index.shtml) (v.1.0.0)
  
 * TE-locate
  * [Perl](http://www.perl.org/get.html) (v.5.14.2)
  * [BWA](http://sourceforge.net/projects/bio-bwa/files/) (v.0.7.2-r351)

 * RetroSeq
  * [Perl](http://www.perl.org/get.html) (v.5.14.2)
  * [BEDTools](https://code.google.com/p/bedtools/downloads/list) (v.2.17.0)
  * [SAMtools](http://sourceforge.net/projects/samtools/files/) (v.0.1.19-44428cd)
  * [BCFTools](https://github.com/samtools/bcftools) (v.0.1.19-44428cd)
  * [Exonerate](http://www.ebi.ac.uk/~guy/exonerate/) (v.2.2.0)
  * [BWA](http://sourceforge.net/projects/bio-bwa/files/) (v.0.7.2-r351)

* TEMP
  * [Perl](http://www.perl.org/get.html) (v.5.14.2)
  * [BioPerl](http://www.bioperl.org/wiki/Getting_BioPerl) (v.1.006901)
  * [BWA](http://sourceforge.net/projects/bio-bwa/files/) (v.0.7.2-r351)
  * [SAMtools](http://sourceforge.net/projects/samtools/files/) (v.0.1.19-44428cd)
  * [BEDTools](https://code.google.com/p/bedtools/downloads/list) (v.2.17.0)
  * [twoBitToFa](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa)

How to run 
------

###Installation
To install the software, from the main pipeline folder, first clone the repository:

```
git clone git@github.com:bergmanlab/mcclintock.git
```

Then cd into the project directory and run the script install.sh with no arguments:

```
cd mcclintock
sh install.sh
```

This will download and unpack all of the TE detection pipelines and check that the required dependencies are available in your path. Missing dependencies will be reported and you must install or make sure these are available to run the full pipeline.

###Running on a test dataset
A script is included to run the full pipeline on a test Illumina resequencing dataset from the yeast genome. To run this test script change directory into the folder named test and run the script runttest.sh. 

```
cd test
sh runtest.sh
```

This script will download the UCSC sacCer2 yeast reference genome, an annotation of TEs in the yeast reference genome from [Carr, Bensasson and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0050978), and a pair of fastq files from SRA, then run the full pipeline.

###Running the pipeline
The pipeline is invoked by running the mcclintock.sh script in the main project folder. This script takes the following 6 input files, specified as options, and will run all five TE detection methods:
* -r : A reference genome sequence in fasta format. (Required)
* -c : The consensus sequences of the TEs for the species in fasta format. (Required)
* -g : The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry. (Required)
* -t : A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file. (Required)
* -b : Retain the sorted and indexed BAM file of the paired end data aligned to the reference genome.
* -i : If this option is specified then all sample specific intermediate files will be removed, leaving only the overall results. The default is to leave sample specific intermediate files (may require large amounts of disk space).
* -1 : The absolute path of the first fastq file from a paired end read, this should be named ending _1.fastq. (Required)
* -2 : The absolute path of the second fastq file from a paired end read, this should be named ending _2.fastq. (Required)
* -p : The number of processors to use for parallel stages of the pipeline. (Default is 1)
* -h : Prints this help guide.

Example pipeline run:
```
sh mcclintock.sh -r reference.fasta -c te_consensus.fasta -g te_locations.gff -t te_families.tsv -1 sample_1.fastq -2 sample_2.fastq -p 2
```

Data created during pre-processing will be stored in a folder in the main directory named after the reference genome used with individual sub-directories for samples. 

###Output format
The output of the run scripts is a bed format file with the 4th column containing the name of the TE name and whether it is a novel insertion (new) or a TE shared with the reference (old). The outputs also include a header line for use with the UCSC genome browser. The final results files are located in a results folder saved in the specific sample folder within the directory named after the reference genome.

###Running individual TE detection methods
Each folder contains one of the TE detection methods tested in the review. In addition to the standard software there is also a file named runXXXX.sh. Running this file without arguments will explain to the user what input files should be used to execute the method. These arguments should be supplied after the script name with spaces in between, as follows:
```
sh runXXXX.sh argument1 argument2 argument3 ...
```

