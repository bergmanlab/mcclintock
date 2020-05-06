<p align="center">
    <img src="https://github.com/bergmanlab/mcclintock/blob/master/img/mcclintock.jpg?raw=true" alt="McClintock in action"/>
</p>

Introduction
------
Many methods have been developed to detect transposable element (TE) insertions from whole genome shotgun next-generation sequencing (NGS) data, each of which has different dependencies, run interfaces, and output formats. Here, we have developed a meta-pipeline to install and run six methods for detecting TE insertions in NGS data, which generates output in the UCSC Browser extensible data (BED) format. A detailed description of the McClintock pipeline and evaluation of McClintock component methods on the yeast genome can be found in [Nelson, Linheiro and Bergman (2017) *G3* 7:2763-2778](http://www.g3journal.org/content/7/8/2763).

The complete pipeline requires a fasta reference genome, a fasta consensus set of TE sequences present in the organism and fastq paired end sequencing reads. Optionally if detailed annotation of the reference TE sequences has been performed, a GFF file with the locations of known TEs present in the reference genome and a tab delimited hierarchy file linking these individual insertions to the consensus they belong to (an example of this file is included in the test directory as sac_cer_te_families.tsv) can be supplied. If only single ended sequencing data are available then this can be supplied as option -1 however only ngs_te_mapper and RelocaTE will run as these are the only methods that handle it.

Software Methods
------
 * [ngs_te_mapper](https://github.com/bergmanlab/ngs_te_mapper) - [Linheiro and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371/journal.pone.0030008)
 * [RelocaTE](https://github.com/srobb1/RelocaTE) - [Robb *et al.* (2013)](http://www.g3journal.org/content/3/6/949.long)
 * [RelocaTE2](https://github.com/stajichlab/RelocaTE2) - [Chen *et al.* (2017)](https://peerj.com/articles/2942/)
 * [TEMP](https://github.com/JialiUMassWengLab/TEMP) - [Zhuang *et al.* (2014)](http://nar.oxfordjournals.org/content/42/11/6826.full)
 * [RetroSeq](https://github.com/tk2/RetroSeq) - [Keane *et al.* (2012)](http://bioinformatics.oxfordjournals.org/content/29/3/389.long)
 * [PoPoolationTE](https://sourceforge.net/projects/popoolationte/) - [Kofler *et al.* (2012)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002487;jsessionid=2CFC9BF7DEF785D90070915204B5F846)
 * [TE-locate](https://sourceforge.net/projects/te-locate/) - [Platzer *et al.* (2012)](http://www.mdpi.com/2079-7737/1/2/395)

Software Dependencies
------
This system was designed to run on linux operating systems. Installation of software dependencies requires conda (see installation)


Installation
-------
Conda is required to install software dependencies. You can install conda via miniconda.

1. Install python3 miniconda (miniconda is a lightweight installer for the conda package manager).
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
```
2. Update conda
```
conda update conda
```
3. Install and activate mcclintock main environment. This will create an environment containing the dependencies require to run the main script (Python3, SnakeMake)
```
git clone git@github.com:bergmanlab/mcclintock.git
cd mcclintock
git checkout refactor ## OMIT WHEN REFACTOR BRANCH IS MERGED WITH MASTER
conda env create -f install/envs/mcclintock.yml --name mcclintock
conda activate mcclintock
```
**NOTE: the `mcclintock` conda environment must always be activate when executing the mcclintock.py script.**

4. Create component method conda environments and install their depenencies. Each software method has their own set of dependencies. To avoid dependency clashes, McClintock creates a separate conda environment for each method and installs their dependencies separately.
```
python3 mcclintock.py --install
```

How to run 
------
### Running on a test dataset
Some test data is provided in the `test/` directory, though the larger files (reads, reference genome) must be downloaded using the `test/download_test_data.sh` script.

```
sh test/download_test_data.sh
```

The test data provided is a UCSC sacCer2 yeast reference genome, an annotation of TEs in the yeast reference genome from [Carr, Bensasson and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0050978), and a pair of fastq files from SRA.

McClintock can then be run with the test data using the following command (make sure mcclintock conda environment is active: `conda active mcclintock`) 

```
python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -g test/reference_TE_locations.gff \
    -t test/sac_cer_te_families.tsv \
    -1 test/SRR800842_1.fastq \
    -2 test/SRR800842_2.fastq \
    -p 4 \
    -o /path/to/output/directory
```
  * change `/path/to/output/directory` to a real path where you desire the mcclintock output to be created.
  * you can also increase `-p 4` to a higher number if you have more CPU threads available.


### Running the pipeline
The pipeline is invoked by running the `mcclintock.py` script in the main project directory. The options available for this script are as follows:
```
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        A reference genome sequence in fasta format
  -c CONSENSUS, --consensus CONSENSUS
                        The consensus sequences of the TEs for the species in
                        fasta format
  -1 FIRST, --first FIRST
                        The path of the first fastq file from paired end read
                        sequencing or the fastq file from single read
                        sequencing
  -2 SECOND, --second SECOND
                        The path of the second fastq file from a paired end
                        read sequencing
  -p PROC, --proc PROC  The number of processors to use for parallel stages of
                        the pipeline [default = 1]
  -M MEM, --mem MEM     The amount of memory available for the pipeline in GB.
                        [default = 4]
  -o OUT, --out OUT     An output folder for the run. [default = '.']
  -m METHODS, --methods METHODS
                        A comma-delimited list containing the software you
                        want the pipeline to use for analysis. e.g. '-m
                        relocate,TEMP,ngs_te_mapper' will launch only those
                        three methods
  -g LOCATIONS, --locations LOCATIONS
                        The locations of known TEs in the reference genome in
                        GFF 3 format. This must include a unique ID attribute
                        for every entry
  -t TAXONOMY, --taxonomy TAXONOMY
                        A tab delimited file with one entry per ID in the GFF
                        file and two columns: the first containing the ID and
                        the second containing the TE family it belongs to. The
                        family should correspond to the names of the sequences
                        in the consensus fasta file
  -s COVERAGE_FASTA, --coverage_fasta COVERAGE_FASTA
                        A fasta file that will be used for TE-based coverage
                        analysis, if not supplied then the consensus sequences
                        of the TEs will be used for the analysis
  -T, --comments        If this option is specified then fastq comments (e.g.
                        barcode) will be incorporated to SAM output. Warning:
                        do not use this option if the input fastq files do not
                        have comments
  -C, --include_consensus_te
                        This option will include the consensus TE sequences as
                        extra chromosomes in the reference file (useful if the
                        organism is known to have TEs that are not present in
                        the reference strain)
  -R, --include_reference_te
                        This option will include the reference TE sequences as
                        extra chromosomes in the reference file
  --clean               This option will make sure mcclintock runs from
                        scratch and doesn't reuse files already created
  --install             This option will install the dependencies of
                        mcclintock
  --debug               This option will allow snakemake to print progress to
                        stdout

```

Example pipeline run:
```
python3 mcclintock.py \
    -r sacCer2.fasta \
    -c sac_cer_TE_seqs.fasta \
    -g reference_TE_locations.gff \
    -t sac_cer_te_families.tsv \
    -1 SRR800842_1.fastq \
    -2 SRR800842_2.fastq \
    -p 20 \
    -m temp,retroseq,trimgalore \
    -o .
```
* The above example will run TEMP, RetroSeq, and Trim_Galore on the data specified. There are other component methods available in the mcclintock pipeline
```
ngs_te_mapper,relocate,relocate2,temp,retroseq,popoolationte,te-locate,coverage,trimgalore
```
* `trimgalore` and `coverage` are not TE detection methods, but rather optional supplementary components that can be used if desired. 
  * `trimgalore` (https://github.com/FelixKrueger/TrimGalore) performs QC of the fastq files and trims adapters. 
  * `coverage` produces read coverage plots of each TE.
# UPDATE BELOW

### Output
Output files for a full McClintock run or individual components are located in a directory named with the following path structure `/reference_genome/sample_name/results/`. All directories referred to below are contained within this parent directory.

If FastQC was present on the system, then output of FastQC will be stored in the directory `/summary/fastqc_analysis`. The `/summary` directory also includes a file with the output of bamstats, a file with the estimated median insert size of the library for paired-end samples, a file with the summary report for the McClintock pipeline and a file with the family and method based TE detection results.

The summary report for the McClintock pipeline `summary_report.txt` will include version and arguments used for the pipeline, number of reads, median insert size, average genome coverage and total number of TEs detected for each method (including all, reference and non-reference).

The file structure for `summary` will look like this:

File | Description
-- | --
/summary/median_insertsize | Median insert size of the library
/summary/bwamem_bamstats.txt | bamstats report
/summary/fastqc_analysis | FastQC analysis output directory
/summary/summary_report.txt | McClintock summary report
/summary/te_summary_table.csv | Family and method based TE detection summary
/summary/te_coverage/te_depth.csv | Normalized average depth for every TE in the TE library (only generated when `-d` or `-D` is specified)

If TE coverage analysis module was enabled using `-d` or `-D` option, then output of the module will be stored in `/te_coverage/te_depth.csv`, containing normalized average depth for every TE in the TE library. Normalized per-base coverage for every TE family will be plotted and saved in `/te_coverage/te_cov_plot` folder.

Original result files for non-reference and (when available) reference predictions produced by each method are stored in the directory `/originalmethodresults`, containing one sub-directory for each method. Original non-reference and (when available) reference files are then merged and filtered to create 0-based .bed6 format files. For TEMP, Retroseq, and PoPoolationTE, original files are converted into \*raw.bed files that contain all unfiltered predictions for non-reference predictions plus reference TE predictions if provided by the method. Score (Retroseq, ≥6), read support (TEMP and PopoolationTE, reads found for both ends)  and sample frequency (TEMP, >10%) are used to filter raw prediction to create \*redundant.bed files. For RelocaTE, original files from individual families are merged and used to create \*redundant.bed files. For TE-locate, original files are used to create \*redundant.bed directly. \*redundant.bed are filtered to remove predictions that have the identical coordinates for different TE families, with the prediction having the greatest read support (RelocaTE, TEMP, PopoolationTE, TE-locate) or call status (Retroseq) being retained. Non-identical overlapping predictions are retained. The final output file of the run scripts for individual methods after merging reference and non-reference predictions and filtering is \*nonredundant.bed. For ngs_te_mapper, no filtering or merging is needed, and thus \*nonredundant.bed files are created directly.

Individual methods output predictions in different annotation frameworks, which are standardized in \*nonredundant.bed as follows: 
 * For ngs_te_mapper, non-reference TEs are annotated as 0-based TSDs. Reference TEs are annotated as 0-based intervals which are inferred from the data, not from the reference TE annotation. Strand information is present for both non-reference and reference TEs. No filtering or coordinate conversions are needed to create 0-based predictions in \*nonredundant.bed files. 
 * For RelocaTE, non-reference TEs are annotated as 1-based TSDs and reference TEs are annotated on the same 1-based intervals provided in the reference annotation. Strand information is present for both non-reference and reference TEs. Start coordinates are transformed by -1 to create 0-based predictions in \*nonredundant.bed files. 
 * For TEMP, non-reference TEs with or without split-read support are annotated as 1-based intervals in the Start and End columns on the basis of closest supporting read-pairs. Non-reference TEs with split-read support at both ends are annotated as 1-based TSDs in the Junction1 and Junction2 columns. Strand information is present for non-reference TEs. Non-reference predictions with read support on both ends ("1p1") and a sample frequency of >10% are retained in the \*nonredundant.bed files. For non-reference predictions with split-read support on both ends, 0-based coordinates of TSDs were used in the \*nonredundant.bed files. For all other non-reference predictions, 0-based coordinates of read-pair intervals were used in the *nonredundant.bed files. TEMP does not directly detect whether reference TEs are present in the sample, however the results of the TEMP absence module can be used to infer complementary information about the presence of reference annotated TEs. TEMP reference predictions are labelled with "nonab", representing a TE inferred from non-absence to distinguish them from reference TE predictions based on direct evidence.
 * For Retroseq, non-reference TEs are annotated as 1-based intervals in the POS column of the VCF file and two consecutive coordinates in the INFO field. No predictions are made for reference TEs. Strand information is not provided. Insertion with a call status of ≥6 are retained, and the two consecutive coordinates in the INFO field are used to represent a single base insertion on 0-based coordinates in the *redundant.bed and \*nonredundant.bed files.
 * For PoPoolationTE, non-reference and reference TEs are annotated as 1-based intervals on either end of the predicted insertion, and also as a midpoint between the inner coordinates of the two terminal spans (which can lead to half-base midpoint coordinates). Intervals of reference TEs are inferred from the data using locations of reference TE annotations.  Only predictions with read support on both ends ("FR") were retained, and the entire interval between the inner coordinates of the of the two terminal spans (not the midpoint) was converted to 0-based coordinates and used in the \*redundant.bed and \*nonredundant.bed files.
 * For TE-locate, non-reference ("new") and reference ("old") TEs are annotated as 1-based positions plus the length of the inserted TE. Strand information is present, where available, for non-reference TEs. For non-reference TEs the single 1-based coordinate is decreased by one for the start of the insertion and used directly for the end on 0-based coordinates in the \*redundant.bed and \*nonredundant.bed files. For reference TEs the single 1-based coordinate is decreased by one for the start of the insertion and increased by the length of the inserted TE for the end on 0-based coordinates in the \*redundant.bed and \*nonredundant.bed files. 
 
The 4th column of \*nonredundant.bed files contains a string with the name of the TE family, whether it is a non-reference or reference insertion, the allele frequency (if available), the sample name, the method name, the evidence type (sr = split-read, rp = read pair), and a numerical index of the prediction in the \*raw.bed or \*redundant.bed file. As noted above, reference insertions for TEMP are labeled as "nonab" to distinguish them from reference TE predictions based on direct evidence. The final \*nonredundant.bed output file also includes a header line for use with the UCSC genome browser. 

License
------
Copyright 2014-2018 Michael G. Nelson, Shunhua Han, and Casey M. Bergman

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
