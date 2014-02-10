Introduction
------
Many methods have been developed to detect TE insertions from whole genome shotgun sequencing data. In this pipeline we have collected five available methods. Most of them have different inputs and output formats so a pipeline has been created that standardises the input and output requirements and returns results from all five methods.

Software Methods
------
 * [ngs_te_mapper](https://github.com/bergmanlab/ngs_te_mapper "Click to go to download location") - [Linheiro and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371/journal.pone.0030008 "Click to go to paper location")
 * [PoPoolationTE](https://code.google.com/p/popoolationte/ "Click to go to download location") - [Kofler *et al.* (2012)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002487;jsessionid=2CFC9BF7DEF785D90070915204B5F846 "Click to go to paper location")
 * [RelocaTE](https://github.com/srobb1/RelocaTE "Click to go to download location") - [Robb *et al.* (2013)](http://www.g3journal.org/content/3/6/949.long "Click to go to paper location")
 * [TE-locate](http://zendto.gmi.oeaw.ac.at/pickup.php?claimID= Y3tZVfN5xipYyBDN&claimPasscode=NArXMbTjmkorWjSM&emailAddr=te_locate%40gmx.at "Click to go to download location") - [Platzer *et al.* (2012)](http://www.mdpi.com/2079-7737/1/2/395 "Click to go to paper location")
 * [RetroSeq](https://github.com/tk2/RetroSeq "Click to go to download location") - [Keane *et al.* (2012)](http://bioinformatics.oxfordjournals.org/content/29/3/389.long "Click to go to paper location")

Software Dependencies
------
All of the software systems must be run on a unix based system with the software dependencies listed per method below. The versions used to run this pipeline are indicated in parentheses and no guarantee is made that it will function using alternate versions.

 * ngs_te_mapper
  * [R](http://cran.r-project.org/) (v.3.0.2)
  * [blat](http://hgwdev.cse.ucsc.edu/~kent/src/blatSrc.zip) (v.35x1)

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

How to run 
------

###Installation
To install the software, from the main pipeline folder, run the script install.sh with no arguments. This will download and unpack all of the TE detection pipelines and check that the required dependencies are available in your path. Missing dependencies will be report and you must install or make sure these are available to run the full pipeline.

###Running on a test dataset
A script is included to run the full pipeline on a test dataset of the yeast genome. To run this test change directory into the folder named test and run the script runttest.sh. This script will download the yeast genome and a pair of fastq files from SRA, then run the full pipeline, and verify the results.

###Running individual TE detection methods
Each folder contains one of the TE detection methods tested in the review. In addition to the standard software there is also a file named runXXXX.sh Running this file without arguments will explain to the user what input files should be used to execute the method. These arguments should be supplied after the script name with spaces in between, as follows:

    runXXXX.sh argument1 argument2 argument3 ...

###Output format
The output of the run scripts is a bed format file with the 4th column containing the name of the TE name and whether it is a novel insertion (new) or a TE shared with the reference (old). The outputs also include a header line for use with the UCSC genome browser.
