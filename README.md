#  McClintock
### Meta-pipeline to identify transposable element insertions using next generation sequencing data
<p align="center">
    <img src="https://github.com/bergmanlab/mcclintock/blob/master/img/mcclintock.jpg?raw=true" alt="McClintock in action"/>
</p>

## <a name="started"></a> Getting Started
```bash
# INSTALL (Requires Conda)
git clone git@github.com:bergmanlab/mcclintock.git
cd mcclintock
git checkout refactor ## OMIT WHEN REFACTOR BRANCH IS MERGED WITH MASTER
conda env create -f install/envs/mcclintock.yml --name mcclintock
conda activate mcclintock
python3 mcclintock.py --install
python3 test/download_test_data.py

# RUN
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

## Table of Contents
* [Getting Started](#started)
* [Introduction](#intro)
* [Software Methods](#methods)
* [Software Dependencies](#dependency)
* [Installing McClintock](#install)
* [Running McClintock](#run)
* [McClintock Output](#output)


## <a name="intro"></a> Introduction
Many methods have been developed to detect transposable element (TE) insertions from whole genome shotgun next-generation sequencing (NGS) data, each of which has different dependencies, run interfaces, and output formats. Here, we have developed a meta-pipeline to install and run six methods for detecting TE insertions in NGS data, which generates output in the UCSC Browser extensible data (BED) format. A detailed description of the McClintock pipeline and evaluation of McClintock component methods on the yeast genome can be found in [Nelson, Linheiro and Bergman (2017) *G3* 7:2763-2778](http://www.g3journal.org/content/7/8/2763).

The complete pipeline requires a fasta reference genome, a fasta consensus set of TE sequences present in the organism and fastq paired end sequencing reads. Optionally if detailed annotation of the reference TE sequences has been performed, a GFF file with the locations of known TEs present in the reference genome and a tab delimited hierarchy file linking these individual insertions to the consensus they belong to (an example of this file is included in the test directory as sac_cer_te_families.tsv) can be supplied. If only single ended sequencing data are available then this can be supplied as option -1 however only ngs_te_mapper and RelocaTE will run as these are the only methods that handle it.

## <a name="methods"></a> Software Methods
 * [ngs_te_mapper](https://github.com/bergmanlab/ngs_te_mapper) - [Linheiro and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371/journal.pone.0030008)
 * [RelocaTE](https://github.com/srobb1/RelocaTE) - [Robb *et al.* (2013)](http://www.g3journal.org/content/3/6/949.long)
 * [RelocaTE2](https://github.com/stajichlab/RelocaTE2) - [Chen *et al.* (2017)](https://peerj.com/articles/2942/)
 * [TEMP](https://github.com/JialiUMassWengLab/TEMP) - [Zhuang *et al.* (2014)](http://nar.oxfordjournals.org/content/42/11/6826.full)
 * [RetroSeq](https://github.com/tk2/RetroSeq) - [Keane *et al.* (2012)](http://bioinformatics.oxfordjournals.org/content/29/3/389.long)
 * [PoPoolationTE](https://sourceforge.net/projects/popoolationte/) - [Kofler *et al.* (2012)](http://www.plosgenetics.org/article/info%3Adoi%2F10.1371%2Fjournal.pgen.1002487;jsessionid=2CFC9BF7DEF785D90070915204B5F846)
 * [TE-locate](https://sourceforge.net/projects/te-locate/) - [Platzer *et al.* (2012)](http://www.mdpi.com/2079-7737/1/2/395)

## <a name="dependency"></a> Software Dependencies
This system was designed to run on linux operating systems. Installation of software dependencies is automated by Conda, thus Conda is required to install McClintock. Conda can be installed via the [Miniconda installer](https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh).

#### Installing Miniconda (Python 3.X)
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME//miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
conda init
```
* `conda init` requires you to close and open a new terminal before it take effect

#### Update Conda
```bash
conda update conda
```

## <a name="install"></a> Installing McClintock
After installing and updating Conda, McClintock can be installed by: **1.** cloning the repository, **2.** creating the Conda environment, and **3.** running the install script

#### Clone McClintock Repository
```bash
git clone git@github.com:bergmanlab/mcclintock.git
cd mcclintock
git checkout refactor ## OMIT WHEN REFACTOR BRANCH IS MERGED WITH MASTER
```

#### Create McClintock Conda Environment
```bash
conda env create -f install/envs/mcclintock.yml --name mcclintock
```
* This installs the base dependencies (`Snakemake`, `Python3`, `BioPython`) needed to run the main mcclintock script into the mcclintock Conda environment

#### Activate McClintock Conda Environment
```bash
conda activate mcclintock
```
* This adds the dependencies installed in the mclintock conda environment to the environment `PATH` so that they can be used by the McClintock scripts.
* **This environment must <ins>always</ins> be activated prior to running any of the McClintock scripts**
* *NOTE: Sometimes activating conda environments does not work via `conda activate myenv` when run through a script submitted to a queuing system, this can be fixed by activating the environment in the script as shown below*
```bash
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate mcclintock
```
* For more on Conda: see the [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/index.html)

#### Install McClintock Component Methods
```
python3 mcclintock.py --install
```
* This command installs each of the TE insertion detection tools and installs a conda environment for each method.

## <a name="run"></a> Running McClintock
Some test data is provided in the `test/` directory, though the fastQ files must be downloaded using the `test/download_test_data.py` script.
```
python test/download_test_data.py
```
* The test data provided is a UCSC sacCer2 yeast reference genome, an annotation of TEs in the yeast reference genome from [Carr, Bensasson and Bergman (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0050978), and a pair of fastq files from SRA.
  
#### Run McClintock with test data
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

#### McClintock Parameters
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

* Available methods to use with `-m/--methods`:
  * `trimgalore` : Runs [Trim Galore](https://github.com/FelixKrueger/TrimGalore) to QC the fastq file(s) and trim the adaptors prior to running the component methods
  * `coverage` : Creates coverage plots for each TE in the fasta provided by `-c/--consensus` or `-s/coverage_fasta` if provided
  * `ngs_te_mapper` : Runs the [ngs_te_mapper](https://github.com/bergmanlab/ngs_te_mapper) component method
  * `relocate` : Runs the [RelocaTE](https://github.com/srobb1/RelocaTE) component method
  * `relocate2` : Runs the [RelocaTE2](https://github.com/stajichlab/RelocaTE2) component method
  * `temp` : Runs the [TEMP](https://github.com/JialiUMassWengLab/TEMP) component method (Paired-End Only)
  * `retroseq` : Runs the [RetroSeq](https://github.com/tk2/RetroSeq) component method (Paired-End Only)
  * `popoolationte` : Runs the [PoPoolation TE](https://sourceforge.net/p/popoolationte/wiki/Main) component method (Paired-End Only)
  * `te-locate` : Runs the [TE-locate](https://sourceforge.net/projects/te-locate) component method (Paired-End Only)

## <a name="output"></a> McClintock Output
The results of the component methods are output to the directory `<output>/results`. 
* Summary files from the run can be located at `<output>/results/summary/`. 
* Each component method has raw output files which can be found at `<output>/results/<method>/unfiltered/`. 
* The results are standardized into a bed format and can be found in `<output>/results/<method>/*.bed` where `<output>/results/<method>/*.nonredundant.bed` has any duplicate predictions removed.
#### Summary files

File | Description
-- | --
`<output>/summary/summary_report.txt` | Summary Report of McClintock run. Contains information on command used, when and where the script was run, and details about the mapped reads
`<output>/summary/te_summary.txt` | a human readable table that shows the number of TE predictions produced from each method (all, reference, non-reference)
`<output>/summary/te_summary.csv` | a comma-delimited table showing TE predictions (all, reference, non-reference) from each method for each TE family 

License
------
Copyright 2014-2018 Michael G. Nelson, Shunhua Han, and Casey M. Bergman

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
