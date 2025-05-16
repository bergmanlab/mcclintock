<p align="center">
    <img src="https://github.com/bergmanlab/mcclintock/blob/4b860571da3aa358b8dbd637831db194ff94ea08/img/mcclintock.jpg?raw=true" alt="McClintock in action"/>
</p>

# **McClintock**: <sub><sup>A meta-pipeline to identify transposable element insertions using short-read whole genome sequencing data</sup></sub>
## <a name="started"></a> Getting Started
```bash
# INSTALL (Requires Conda and Mamba to be installed)
git clone git@github.com:bergmanlab/mcclintock.git
cd mcclintock
mamba env create -f install/envs/mcclintock.yml --name mcclintock
conda activate mcclintock
python3 mcclintock.py --install
python3 test/download_test_data.py

# RUN
python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -g test/reference_TE_locations.gff \
    -t test/sac_cer_te_families.tsv \
    -1 test/SRR800842_1.fastq.gz \
    -2 test/SRR800842_2.fastq.gz \
    -p 4 \
    -o /path/to/output/directory
```

## Table of Contents
* [Getting Started](#started)
* [Introduction](#intro)
* [Installing Conda/Mamba](#conda)
* [Installing McClintock](#install)
* [McClintock Usage](#run)
* [McClintock Input](#input)
* [McClintock Output](#output)
* [Run Examples](#examples)
* [Citation](#citation)
* [License](#license)

## <a name="intro"></a> Introduction
Many methods have been developed to detect transposable element (TE) insertions from short-read whole genome sequencing (WGS) data, each of which has different dependencies, run interfaces, and output formats. McClintock provides a meta-pipeline to reproducibly install, execute, and evaluate multiple TE detectors, and generate output in standardized output formats. A description of the original McClintock 1 pipeline and evaluation of the original six TE detectors on the yeast genome can be found in [Nelson, Linheiro and Bergman (2017) *G3* 7:2763-2778](http://www.g3journal.org/content/7/8/2763). A description of the re-implemented McClintock 2 pipeline, the reproducible simulation system, and evaluation of 12 TE detectors on the yeast genome can be found in [Chen, Basting, Han, Garfinkel and Bergman (2023) *Mobile DNA* 14:8](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-023-00296-4). The set of TE detectors currently included in McClintock 2 are:

 * [ngs_te_mapper](https://github.com/bergmanlab/ngs_te_mapper) - [Linheiro and Bergman (2012)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030008)
 * [ngs_te_mapper2](https://github.com/bergmanlab/ngs_te_mapper2) - [Han *et al.* (2021)](https://academic.oup.com/genetics/article/219/2/iyab113/6321957)
 * [PoPoolationTE](https://sourceforge.net/projects/popoolationte/) - [Kofler *et al.* (2012)](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002487)
 * [PoPoolationTE2](https://sourceforge.net/p/popoolation-te2/wiki/Home) - [Kofler *et al.* (2016)](https://academic.oup.com/mbe/article/33/10/2759/2925581)
 * [RelocaTE](https://github.com/srobb1/RelocaTE) - [Robb *et al.* (2013)](http://www.g3journal.org/content/3/6/949.long)
 * [RelocaTE2](https://github.com/stajichlab/RelocaTE2) - [Chen *et al.* (2017)](https://peerj.com/articles/2942/)
 * [RetroSeq](https://github.com/tk2/RetroSeq) - [Keane *et al.* (2012)](http://bioinformatics.oxfordjournals.org/content/29/3/389.long)
 * [TEBreak](https://github.com/adamewing/tebreak) - [Schauer *et al.* (2018)](https://genome.cshlp.org/content/28/5/639)
 * [TEFLoN](https://github.com/jradrion/TEFLoN) - [Adrion *et al.* (2017)](https://academic.oup.com/gbe/article/9/5/1329/3064433)
 * [TE-locate](https://sourceforge.net/projects/te-locate/) - [Platzer *et al.* (2012)](http://www.mdpi.com/2079-7737/1/2/395)
 * [TEMP](https://github.com/JialiUMassWengLab/TEMP) - [Zhuang *et al.* (2014)](http://nar.oxfordjournals.org/content/42/11/6826.full)
 * [TEMP2](https://github.com/weng-lab/TEMP2) - [Yu *et al.* (2021)](https://academic.oup.com/nar/article/49/8/e44/6123378?login=true)

## <a name="conda"></a> Installing Conda/Mamba via Miniforge
McClintock is written in Python3 leveraging the [SnakeMake](https://snakemake.readthedocs.io/en/stable/) workflow system and is designed to run on linux operating systems. Installation of software dependencies for McClintock and its component methods is automated by [Conda](https://docs.conda.io/en/latest/), thus a working installation of Conda (and it's reimplementation [Mamba](https://mamba.readthedocs.io/en/latest/)) is required to install McClintock. Conda/Mamba can be installed via the [Miniforge installer](https://github.com/conda-forge/miniforge?tab=readme-ov-file#install).

```bash
wget -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3.sh -b -p "${HOME}/conda" 
source "${HOME}/conda/etc/profile.d/conda.sh"
source "${HOME}/conda/etc/profile.d/mamba.sh"
conda init
```
* `conda init` requires you to close and open a new terminal before it take effect

## <a name="install"></a> Installing McClintock
After installing and updating Conda/Mamba, McClintock and its component methods can be installed by: **1.** cloning the repository, **2.** creating the Conda environment, and **3.** running the install script.

#### Clone McClintock Repository
```bash
git clone git@github.com:bergmanlab/mcclintock.git
cd mcclintock
```

#### Create McClintock Conda Environment
```bash
mamba env create -f install/envs/mcclintock.yml --name mcclintock
```
* This installs the base dependencies needed to run the main McClintock script (`Snakemake`, `Python3`, `BioPython`) into the `mcclintock` Conda environment.

#### Activate McClintock Conda Environment
```bash
conda activate mcclintock
```
* This adds the dependencies installed in the McClintock conda environment to the environment `PATH` so that they can be used by the McClintock scripts.
* **This environment must <ins>always</ins> be activated prior to running any of the McClintock scripts**
* *NOTE: Sometimes activating conda environments does not work via `conda activate myenv` when run through a script submitted to a queueing system, this can be fixed by activating the environment in the script as shown below*
```bash
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate mcclintock
```
* For more on Conda: see the [Conda User Guide](https://docs.conda.io/projects/conda/en/latest/index.html)

#### Install McClintock Component Methods
* To install all of the component methods and create their associated conda environments, use the following command:
```bash
python3 mcclintock.py --install
```

* If you only want to install specific methods to save space and time, you can specify method(s) using the `-m` flag:
```bash
python3 mcclintock.py --install -m <method1>,<method2>
```

* *NOTE: If you re-run either the full installation or installation of specific methods, the installation script will do a clean installation and remove previously installed components.*

* If you want to install missing methods to an already existing mcclintock installation, you can use the `--resume` flag:
```bash
python3 mcclintock.py --install --resume
```

* *NOTE: If you use the `--resume` flag when installing specific method(s) with `-m`, the installation script will only install the specified method(s) if they haven't previously been installed. Do not use the `--resume` flag if you want to do a clean installation of a specific method.*

## <a name="run"></a> McClintock Usage

Running the complete McClintock pipeline requires a fasta reference genome (option `-r`), a set of TE consensus/canonical sequences present in the organism (option `-c`), and fastq paired-end sequencing reads (options `-1` and `-2`). If only single-end fastq sequencing data are available, then this can be supplied using only option `-1`, however only the TE detectors that handle single-ended data will be run. Optionally, if a detailed annotation of TE sequences in the reference genome has been performed, a GFF file with annotated reference TEs (option `-g`) and a tab-delimited "taxonomy" file linking annotated insertions to their TE family (option `-t`) can be supplied. Example input files are included in the [test](https://github.com/bergmanlab/mcclintock/blob/master/test/) directory. 

```
##########################
##       Required       ##
##########################
  -r, --reference REFERENCE
                        A reference genome sequence in fasta format
  -c, --consensus CONSENSUS
                        The consensus sequences of the TEs for the species in
                        fasta format
  -1, --first FIRST
                        The path of the first fastq file from paired end read
                        sequencing or the fastq file from single read
                        sequencing

##########################
##       Optional       ##
##########################
  -h, --help            show this help message and exit
  -2, --second SECOND
                        The path of the second fastq file from a paired end
                        read sequencing
  -p, --proc PROC       The number of processors to use for parallel stages of
                        the pipeline [default = 1]
  -o, --out OUT         An output folder for the run. [default = '.']
  -m, --methods METHODS
                        A comma-delimited list containing the software you
                        want the pipeline to use for analysis. e.g. '-m
                        relocate,TEMP,ngs_te_mapper' will launch only those
                        three methods. If this option is not set, all methods
                        will be run [options: ngs_te_mapper, ngs_te_mapper2, 
                        relocate, relocate2, temp, temp2, retroseq, 
                        popoolationte, popoolationte2, te-locate, teflon, 
                        coverage, trimgalore, map_reads, tebreak]

  -g, --locations LOCATIONS
                        The locations of known TEs in the reference genome in
                        GFF 3 format. This must include a unique ID attribute
                        for every entry. If this option is not set, a file of 
                        reference TE locations in GFF format will be produced 
                        using RepeatMasker
  -t, --taxonomy TAXONOMY
                        A tab delimited file with one entry per ID in the GFF
                        file and two columns: the first containing the ID and
                        the second containing the TE family it belongs to. The
                        family should correspond to the names of the sequences
                        in the consensus fasta file. If this option is not set, 
                        a file mapping reference TE instances to TE families 
                        in TSV format will be produced using RepeatMasker
  -s, --coverage_fasta COVERAGE_FASTA
                        A fasta file that will be used for TE-based coverage
                        analysis, if not supplied then the consensus sequences
                        of the TEs set by -c/--consensus will be used for the 
                        analysis
  -a, --augment AUGMENT
                        A fasta file of TE sequences that will be included as
                        extra chromosomes in the reference file (useful if the
                        organism is known to have TEs that are not present in
                        the reference strain)
  -k, --keep_intermediate KEEP_INTERMEDIATE
                        This option determines which intermediate files are 
                        preserved after McClintock completes [default: general]
                        [options: minimal, general, methods, <list,of,methods>, 
                        all]
  -s, --sample_name SAMPLE_NAME
                        The sample name to use for output files [default: 
                        fastq1 name]
  -n, --config CONFIG   This option determines which config files to use for 
                        your McClintock run [default: config in McClintock 
                        Repository]
  -v, --vcf VCF         This option determines which format of VCF output will 
                        be created [default: siteonly][options: siteonly,sample]
  --install             This option will install the dependencies of McClintock
  --resume              This option will attempt to use existing intermediate 
                        files from a previous McClintock run
  --debug               This option will allow snakemake to print progress to 
                        stdout
  --serial              This option runs without attempting to optimize thread 
                        usage to run rules concurrently. Each multithread rule 
                        will use the max processors designated by -p/--proc
  --make_annotations    This option will only run the pipeline up to the 
                        creation of the repeat annotations
  --comments            If this option is specified then fastq comments (e.g.
                        barcode) will be incorporated to SAM output. Warning:
                        do not use this option if the input fastq files do not
                        have comments
```

* Available methods to use with `-m/--methods`:
  * `trimgalore` : Runs [Trim Galore](https://github.com/FelixKrueger/TrimGalore) to QC the fastq file(s) and trim the adaptors prior to running the component methods
  * `coverage` : Estimates copy number based on normalized coverage and creates coverage plots for each TE in the fasta provided by `-c/--consensus` or `-s/coverage_fasta` if provided
  * `map_reads` : Maps the reads to the reference genome. This is useful to ensure the BAM alignment file is produced regardless if another method requires it as input
  * `ngs_te_mapper` : Runs the [ngs_te_mapper](https://github.com/bergmanlab/ngs_te_mapper) component method
  * `ngs_te_mapper2`: Runs the [ngs_te_mapper2](https://github.com/bergmanlab/ngs_te_mapper2) component method
  * `relocate` : Runs the [RelocaTE](https://github.com/srobb1/RelocaTE) component method
  * `relocate2` : Runs the [RelocaTE2](https://github.com/stajichlab/RelocaTE2) component method
  * `temp` : Runs the [TEMP](https://github.com/JialiUMassWengLab/TEMP) component method (Paired-End Only)
  * `temp2` : Runs the [TEMP2](https://github.com/weng-lab/TEMP2) component method (Paired-End Only)
  * `retroseq` : Runs the [RetroSeq](https://github.com/tk2/RetroSeq) component method (Paired-End Only)
  * `popoolationte` : Runs the [PoPoolation TE](https://sourceforge.net/p/popoolationte/wiki/Main) component method (Paired-End Only)
  * `popoolationte2` : Runs the [PoPoolation TE2](https://sourceforge.net/p/popoolation-te2/wiki/Home) component method (Paired-End Only)
  * `te-locate` : Runs the [TE-locate](https://sourceforge.net/projects/te-locate) component method (Paired-End Only)
  * `teflon` : Runs the [TEFLoN](https://github.com/jradrion/TEFLoN) component method (Paired-End Only)

## <a name="input"></a> Mcclintock Input Files
#### Warning
  * Feature names (contig IDs, TE IDs, Family IDs) in input files must not contain any invalid symbols to ensure compatibility with all component methods.
  * INVALID_SYMBOLS:`; & ( ) | * ? [ ] ~ { } < ! ^ " , \  $ / + - #`

#### Required
* Reference FASTA (`-r/--reference`)
  * The genome sequence of the reference genome in FASTA format. The reads from the FASTQ file(s) will be mapped to this reference genome to predict TE insertions
  * [example](https://github.com/bergmanlab/mcclintock/blob/master/test/sacCer2.fasta)
* Consensus FASTA (`-c/--consensus`)
  * A FASTA file containing a consensus sequence for each family
  * [example](https://github.com/bergmanlab/mcclintock/blob/master/test/sac_cer_TE_seqs.fasta)
* FASTQ File 1 (`-1/--first`)
  * Either the Read1 FASTQ file from a paired-end sequencing run or the FASTQ file from an unpaired sequencing run
#### Optional
* FASTQ File 2 (`-2/--second`)
  * The Read2 FASTQ file from a paired-end sequencing run. Not required if using unpaired data.
* Locations (`-g/--locations`)
  * A GFF file contianing the locations of the reference TEs.
  * Each annotation should contain an `ID=` attribute that contains a unique identifier that does not match any other annotation.
  * If both the locations GFF and taxonomy TSV are not provided, McClintock will produce them using RepeatMasker with the consensus sequences
  * [example](https://github.com/bergmanlab/mcclintock/blob/master/test/reference_TE_locations.gff)
* Taxonomy (`-t/--taxonomy`)
  * A Tab delimited file that maps the unique reference TE to the family it belongs to.
  * This file should contain two columns, the first corresponding to the reference TE identifier which should match the `ID=` attribute from the locations GFF(`-g`). The second column contains the reference TE's family which should match the name of a sequence in the consensus fasta (`-c`)
  * If both the locations GFF and taxonomy TSV are not provided, McClintock will produce them using RepeatMasker with the consensus sequences
  * [example](https://github.com/bergmanlab/mcclintock/blob/master/test/sac_cer_te_families.tsv)
* Coverage FASTA (`-s/--coverage_fasta`)
  * A fasta file of TE sequences to be used for the coverage analysis.
  * By default, McClintock estimates the coverage and creates coverage plots of the consensus TE sequences (`-c`). This option allows you to use a custom set of TEs for the coverage estimations and plots.
* Augment FASTA (`-a/--augment`)
  * A FASTA file of TE sequences that will be included as extra chromosomes in the reference genome file (`-r`)
  * Some methods leverage the reference TE sequences to find non-reference TEs insertions. The augment FASTA can be used to augment the reference genome with additional TEs that can be used to locate non-reference TE insertions that do not have a representative in the reference genome.

## <a name="output"></a> McClintock Output
The results of McClintock component methods are output to the directory `<output>/<sample>/results`.
* Summary files from the run can be located at `<output>/<sample>/results/summary/`.
* Each component method has raw output files which can be found at `<output>/<sample>/results/<method>/unfiltered/`.
* Raw results are filtered by parameters defined in the `<method>_post.py` postprocessing configuration files for each method, then standardized into BED and VCF formats. Post-processing config files can be found in `/path/to/mcclintock/config/<method>` and can be modified if you want to adjust default filtering parameters. 
* Standardize results in BED format contain non-reference and (when available for a method) reference TE predictions and can be found in `<output>/<sample>/results/<method>/*.bed`, where `<output>/<sample>/results/<method>/*.nonredundant.bed` has any redundant predictions removed.
* Standardized results in VCF format contain only non-reference TE predictions and can be found in  `<output>/<sample>/results/<method>/*_nonredundant_non-reference_*.vcf`. Two types of VCF are supported: (i) VCF files with "site-only" information (`*_nonredundant_non-reference_siteonly.vcf`) and (ii) VCF files that contains a "sample" column (`*_nonredundant_non-reference_sample.vcf`). The position of TE insertion variants in VCF files corresponds to the start postion of predicted intervals in nonredundant BED files. Note that the genotype (GT) field of `*_nonredundant_non-reference_sample.vcf` only indicates the presence/absence of a non-reference TE insertion variant call in the sample, and does not contain information about the ploidy or zygosity of the variant in the sample.

#### HTML Summary Report: `<output>/<sample>/results/summary/`
* McClintock generates an interactive HTML summary report that contains information on how the run was executed, read mapping information, QC information, and a summary of component method predictions. `<output>/<sample>/results/summary/summary.html`
* This page also links to the pages that summarize the predictions from each method: all predictions by method, predictions for each family, predictions for each contig. `<output>/<sample>/results/summary/html/<method>.html`
* The HTML report also summarizes reference and non-reference predictions for all families. `<output>/<sample>/results/summary/html/families.html`
* A page is also generated for each family, which summarizes the coverage for the family consensus sequence and the family-specific predictions from each component method. `<output>/<sample>/results/summary/html/<family>.html`

#### Raw Summary files : `<output>/<sample>/results/summary/`
* `<output>/<sample>/results/summary/data/run/summary_report.txt` : Summary Report of McClintock run. Contains information on the McClintock command used, when and where the script was run, details about the mapped reads, and table that shows the number of TE predictions produced from each method.
* `<output>/<sample>/results/summary/data/run/te_prediction_summary.txt` : A comma-delimited table showing reference and non-reference predictions for each component method
* `<output>/<sample>/results/summary/data/families/family_prediction_summary.txt` : a comma-delimited table showing TE predictions (all, reference, non-reference) from each method for each TE family
* `<output>/<sample>/results/summary/data/coverage/te_depth.txt` : (Only produced if coverage module is run) a comma-delimited table showing normalized depth for each consensus TE or TE provided in coverage fasta.
* All tables and plots contain a link to the raw data so that users can manually filter or visualize it with other programs.

#### TrimGalore : `<output>/<sample>/results/trimgalore/`
* `<fastq>_trimming_report.txt` : Information on parameters used and statistics related to adapter trimming with cutadapt. Provides an overview of sequences removed via the adapter trimming process.
* `<fastq>_fastqc.html` : FastQC report of the trimmed fastq files. Provides information on the results of steps performed by FastQC to assess the quality of the trimmed reads.
* `<fastq>_fastqc.zip` : FastQC graph images and plain-text summary statistics compressed into a single `.zip` file

#### Coverage : `<output>/<sample>/results/coverage/`
* `plots/*.png` : Coverage plots showing the normalized read coverage across each TE either from the consensus fasta (`-c`) or the coverage fasta (`-s`) if provided. Coverage of uniquely mapping reads (MAPQ > 0) is in dark gray, while coverage of all reads (MAPQ >= 0) is in light gray. Raw coverage at each postion in a TE is normalized to the average mapping depth at unique regions of the hard-masked reference genome. The average normalized coverage is shown as a black line, and is estimated from the central region of each TE omitting regions at the 5' and 3' ends equal to the average read length to prevent biases due to  mapping at TE edges.
* `te-depth-files/*.allQ.cov` : Raw read coverage at each position in a TE sequence. (Output of `samtools depth`)
* `te-depth-files/*.highQ.cov` : coverage of mapped reads with MAPQ > 0 at each position, omitting multi-mapped reads.

#### ngs_te_mapper : `<output>/<sample>/results/ngs_te_mapper/`
* `unfiltered/<sample>_insertions.bed` : BED file containing raw 0-based intervals corresponding to TSDs for non-reference predictions and 0-based intervals corresponding to the reference TEs. Reference TE intervals are inferred from the data, not from the reference TE annotations. Strand information is present for both non-reference and reference TEs.
* `<sample>_ngs_te_mapper_nonredundant.bed` : BED file containing 0-based intervals corresponding to TSDs for non-reference predictions and 0-based intervals corresponding to the reference TEs. This file contains the same predictions from `unfiltered/<sample>_insertions.bed` with the BED line name adjusted to match the standard McClintock naming convention. By default, no filtering is performed on the raw `ngs_te_mapper` predictions aside from removing redundant predictions. However, the config file: (`/path/to/mcclintock/config/ngs_te_mapper/ngs_te_mapper_post.py`) can be modified to increase the minimum read support threshold if desired.

#### ngs_te_mapper2 : `<output>/<sample>/results/ngs_te_mapper2/`
* `unfiltered/<sample>.nonref.bed`: BED file containing raw 0-based intervals corresponding to the 5' and 3' breakpoints for non-reference predictions.
* `unfiltered/<sample>.ref.bed`: BED file containing raw 0-based intervals corresponding to the reference TE annotations predicted by ngs_te_mapper2
* `<sample>_ngs_te_mapper2_nonredundant.bed`: BED file containing predictions from `unfiltered/<sample>.nonref.bed` and `unfiltered/<sample>.ref.bed` with BED line names matching the standard McClintock naming convention.

#### PoPoolationTE : `<output>/<sample>/results/popoolationTE/`
* `unfiltered/te-poly-filtered.txt` : Tab-delimited table with non-reference and reference TE predictions and support values. Predictions are annotated as 1-based intervals on either end of the predicted insertion, and also as a midpoint between the inner coordinates of the two terminal spans (which can lead to half-base midpoint coordinates)
* `<sample>_popoolationte_nonredundant.bed` : BED file containing only TE predictions with read support on both ends ("FR") and with percent read support >0.1 for both ends were retained in this file. The entire interval between the inner coordinates of the of the two terminal spans (not the midpoint) was converted to 0-based coordinates. Filtering parameters can be modified in the config file: (`/path/to/mcclintock/config/popoolationte/popoolationte_post.py`)

#### PoPoolationTE2 : `<output>/<sample>/results/popoolationTE2/`
* `unfiltered/teinsertions.txt` : Tab-delimited table with TE predictions and TE frequency values (ratio of physical coverage supporting a TE insertion to the total physical coverage). PoPoolationTE2 does not indicate which predictions are reference and non-reference TEs. Also, only a single position is reported for each prediction, so the TSD is not predicted. Predictions may only have support from one side of the junction ("F" or "R") or both sides ("FR"). Prediction coordinates are 1-based.
* `<sample>_popoolationte2_nonredundant.bed` : BED file containing all of the predictions from `unfiltered/teinsertions.txt` that have support on both ends ("FR") and have a frequency >0.1. The filtering criteria can be modified in the config file: (`/path/to/mcclintock/config/popoolationte2/popoolationte2_post.py`). If predictions overlap a TE in the reference genome, that reference TE is reported in this file using the positions of the reference TE annotation (not the position reported by PoPoolationTE2). If the prediction does not overlap a reference TE, it is designated a non-reference TE insertion `_non-reference_`. The coordinates for all predictions are adjusted to be 0-based.

#### RelocaTE : `<output>/<sample>/results/relocaTE/`
* `unfiltered/combined.gff` : GFF containing 1-based TSDs for non-reference predictions and 1-based intervals for reference TEs. The reference intervals are based on the reference TE annotations.
* `<sample>_relocate_nonredundant.bed` : BED file containing predictions from `unfiltered/combined.gff` converted into 0-based intervals with BED line names matching the standard McClintock naming convention. By default, no filtering is performed on the raw predictions aside from removing redundant predictions. However, the config file: (`/path/to/mcclintock/config/relocate/relocate_post.py`) can be modified to increase the minimum left and right prediction support thresholds for both reference and non-reference predictions.

#### RelocaTE2 : `<output>/<sample>/results/relocaTE2/`
* `unfiltered/repeat/results/ALL.all_ref_insert.gff` : GFF file containing reference TE predictions with 1-based coordinates. The final column also contains read counts supporting the junction (split-read) and read counts supporting the insertion (read pair).
* `unfiltered/repeat/results/ALL.all_nonref_insert.gff` : GFF file containing non-reference TE predictions with 1-based coordinates. The final column also contains read counts supporting the junction (split-read) and read counts supporting the insertion (read pair).
* `<sample>_relocate2_nonredundant.bed` : BED file containing all reference and non-reference predictions from `ALL.all_ref_insert.gff` and `ALL.all_nonref_insert.gff`. Coordinates are adjusted to be 0-based. By default, no filtering is performed on split-read and split-pair evidence. However, the config file: (`/path/to/mcclintock/config/ngs_te_mapper/ngs_te_mapper_post.py`) can be modified to increase the default threshold for these values.

#### RetroSeq : `<output>/<sample>/results/retroseq/`
* `unfiltered/<sample>.call` : VCF file containing non-reference TE predictions. Non-reference TEs are annotated as 1-based intervals in the POS column and two consecutive coordinates in the INFO field. No predictions are made for reference TEs. Strand information is not provided.
* `<sample>_retroseq_nonredundant.bed` : BED file containing non-reference TE predictions from `unfiltered/<sample>.call` with a [Breakpoint confidence threshold](https://github.com/tk2/RetroSeq/wiki/RetroSeq-Tutorial#interpreting-the-output) of >6 are retained in this file. This filtering threshold can be changed by modifying the config file: (`/path/to/mcclintock/config/retroseq/retroseq_post.py`). The position interval reported in the `INFO` column of `unfiltered/<sample>.call` VCF file are converted to 0-based positions and used as the start and end positions in the BED lines in this file.

#### TEbreak : `<output>/<sample>/results/tebreak/`
* `unfiltered/<sample>.sorted.tebreak.table.txt` : Tab-delimited table containing non-reference TE predictions with 0-based coordinates. No predictions are made for reference TEs. Strand information is provided.
* `<sample>_tebreak_nonredundant.bed` : BED file containing non-reference TE predictions from `unfiltered/<sample>.sorted.tebreak.table.txt` with no additional filtering. This filtering threshold can be changed by modifying the config file: (`/path/to/mcclintock/config/retroseq/tebreak_post.py`).

#### TEMP : `<output>/<sample>/results/TEMP/`
* `unfiltered/<sample>.absence.refined.bp.summary` : Tab-delimited table containing reference TEs that are predicted to be absent from the short read data. Position intervals are 1-based.
* `unfiltered/<sample>.insertion.refined.bp.summary` : Tab-delimited table containing non-reference TE predictions. Position intervals are 1-based.
* `<sample>_temp_nonredundant.bed` : BED file containing all reference TEs not reported as absent by TEMP in the `unfiltered/<sample>.absence.refined.bp.summary` file. Also contains non-reference TE predictions `unfiltered/<sample>.insertion.refined.bp.summary` formatted as a BED line using the McClintock naming convention. Positions for both reference and non-reference predictions are 0-based. Non-reference predictions from `unfiltered/<sample>.insertion.refined.bp.summary` are only added to this file if the prediction has read support on both ends ("1p1") and has a sample frequency of > 10%. These filtering restrictions can be modified in the config file: (`/path/to/mcclintock/config/TEMP/temp_post.py`). Non-reference TEs with split-read support at both ends are marked in the BED line name with "_sr_" and the `Junction1` and `Junction2` columns from `unfiltered/<sample>.insertion.refined.bp.summary` are used as the start and end positions of the TSD in this file (converted to 0-based positions). If the non-reference TE prediction does not have split-read support on both ends of the insertions, the designation "_rp_" is added to the BED line name and the `Start` and `End` columns from `unfiltered/<sample>.insertion.refined.bp.summary` are used as the start and end positions of the TSD in this file (converted to 0-based). Note: TEMP reference insertions are labeled `nonab` in the BED line name since they are inferred by no evidence of absence, in contrast to reference insertions detected by other components that are inferred from direct evidence of their presence.

#### TEMP2 : `<output>/<sample>/results/temp2/`
* `unfiltered/<sample>.absence.refined.bp.summary` : Tab-delimited table containing reference TEs that are predicted to be absent from the short read data. Position intervals are 1-based.
* `unfiltered/<sample>.insertion.bed`: BED file containing 0-based coordinates of the non-reference predictions.
* `<sample>_temp_nonredundant.bed`: BED file containing all reference TEs not reported in the `unfiltered/<sample>.absence.refined.bp.summary`. Also contains non-reference TE predictions from `unfiltered/<sample>.insertion.bed`. Non-reference predictions from `unfiltered/<sample>.insertion.bed` are only added to this file if the prediction has read support on both ends ("1p1") and has a sample frequency of > 10%. These filtering restrictions can be modified in the config file: (`/path/to/mcclintock/config/temp2/temp2_post.py`).

#### TE-Locate : `<output>/<sample>/results/te-locate/`
* `unfiltered/te-locate-raw.info` : A tab-delimited table containing reference ("old") and non-reference ("new") predictions using 1-based positions. TSD intervals are not predicted for non-reference TEs, instead a single position is reported.
* `<sample>_telocate_nonredundant.bed` : BED file containing all reference and non-reference predictions from `unfiltered/te-locate-raw.info`. Coordinates for both reference and non-reference TE predictions are converted to a 0-based interval. The reference TE end position is extended by the `len` column in `unfiltered/te-locate-raw.info`. Non-reference TE predictions are a single position as TE-Locate does not predict the TSD size.

#### TEFLoN : `<output>/<sample>/results/teflon/`
* `unfiltered/genotypes/sample.genotypes.txt` : A tab-delimited table containing all of the breakpoints and support information for insertion predictions. Predictions are treated as reference predictions if they contain a TE ID in column 7.
```bash
# from: https://github.com/jradrion/TEFLoN
C1: chromosome
C2: 5' breakpoint estimate ("-" if estimate not available)
C3: 3' breakpoint estimate ("-" if estimate not available)
C4: search level id (Usually TE family)
C5: cluster level id (Usually TE order or class)
C6: strand ("." if strand could not be detected)
C7: reference TE ID ("-" if novel insertion)
C8: 5' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
C9: 3' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")
C10: read count for "presence reads"
C11: read count for "absence reads"
C12: read count for "ambiguous reads"
C13: genotype for every TE (allele frequency for pooled data, present/absent for haploid, present/absent/heterozygous for diploid) #Note: haploid/diploid caller is under construction, use "pooled" for presence/absence read counts
C14: numbered identifier for each TE in the population
```
* `<sample>_teflon_nonredundant.bed` : BED file containing all reference and non-reference predictions from `unfiltered/genotypes/sample.genotypes.txt`. Reference predictions use the coordinates for the TE with the reference ID from column 7. By default, only non-reference predictions with both breakpoints (C2 and C3) are kept in this file. Non-reference predictions must also have at least 3 presence reads (C10) and an allele frequency greater than `0.1` (C13). These filtering restrictions can be changed by modifying the TEFLoN config file: `/path/to/mcclintock/config/teflon/teflon_post.py`

## <a name="examples"></a> Run Examples
#### Running McClintock with test data
This repository also provides test data to ensure your McClintock installation is working. Test data can be found in the `test/` directory which includes a yeast reference genome (UCSC sacCer2) and an annotation of TEs in this version of the yeast genome from [Carr et al. (2012)](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0050978). A pair of fastq files can be downloaded from SRA using the `test/download_test_data.py` script:
```
python test/download_test_data.py
```
* Once the fastq files have been downloaded, Mcclintock can be run on the test data as follows:
```
python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -g test/reference_TE_locations.gff \
    -t test/sac_cer_te_families.tsv \
    -1 test/SRR800842_1.fastq.gz \
    -2 test/SRR800842_2.fastq.gz \
    -p 4 \
    -o /path/to/output/directory
```
* Change `/path/to/output/directory` to a real path where you desire the McClintock output to be created.
* You can also increase `-p 4` to a higher number if you have more CPU threads available.
* A working installation of McClintock applied to the test data should yield the following results in `/path/to/output/directory/SRR800842_1/results/summary/data/run/summary_report.txt`:

```
----------------------------------
MAPPED READ INFORMATION
----------------------------------
read1 sequence length:  94
read2 sequence length:  94
read1 reads:            18547818
read2 reads:            18558408
median insert size:     302
avg genome coverage:    268.185
----------------------------------

-----------------------------------------------------
METHOD          ALL       REFERENCE    NON-REFERENCE 
-----------------------------------------------------
ngs_te_mapper   35        21           14            
ngs_te_mapper2  87        49           38            
relocate        80        63           17            
relocate2       139       41           98            
temp            365       311          54            
temp2           367       311          56            
retroseq        58        0            58            
popoolationte   141       130          11            
popoolationte2  186       164          22            
te-locate       713       164          549           
teflon          414       390          24            
tebreak         60        0            60            
-----------------------------------------------------
```
* NOTE: `popoolationte` and `popoolationte2` exhibit run-to-run variation so numbers for these methods will differ slightly on replicate runs of the test data.

#### Running McClintock with specific component methods
* By default, McClintock runs all components (all 12 TE detection methods plus the coverage module using the output of the trimgalore method).
* If you only want to run a specific component method, you can use the `-m` flag to specify which method to run:
```
python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -g test/reference_TE_locations.gff \
    -t test/sac_cer_te_families.tsv \
    -1 test/SRR800842_1.fastq.gz \
    -2 test/SRR800842_2.fastq.gz \
    -p 4 \
    -m temp \
    -o /path/to/output/directory
```
* You can also specify an arbitrary set of multiple component methods to run using a comma-separated list of the methods after the `-m` flag:
```
python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -g test/reference_TE_locations.gff \
    -t test/sac_cer_te_families.tsv \
    -1 test/SRR800842_1.fastq.gz \
    -2 test/SRR800842_2.fastq.gz \
    -p 4 \
    -m trimgalore,temp,ngs_te_mapper,retroseq \
    -o /path/to/output/directory
```
* Note: if the `-m` option is set, you must specify the `trimgalore` method explicitly if you want the other component methods to use trimmed reads as input.

#### Running McClintock with multiple samples using same reference genome
* When running McClintock on multiple samples that use the same reference genome and consensus TEs, it is advised to generate the reference TE annotation GFF and TE Taxonomy TSV files as a pre-processing step. Otherwise, these files will be created by McClintock for each sample, which can lead to increased run time and disk usage.
* To create the reference TE annotation GFF and TE Taxonomy TSV files, you can run McClintock with the `--make_annotations` flag, which will use RepeatMasker to produce only these files, then exit:
```
python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -p 4 \
    -o <output> \
    --make_annotations
```
* The locations of the reference TE annotation GFF and TE Taxonomy TSV files generated using `--make_annotations` are as follows:
  * Reference TE locations GFF: `<output>/<reference_name>/reference_te_locations/unaugmented_inrefTEs.gff`
  * TE Taxonomy TSV:  `<output>/<reference_name>/te_taxonomy/unaugmented_taxonomy.tsv`
* You can then use the `--resume` flag for future runs with the same reference genome and output directory without having to redundantly generate these files for each run:

```bash
python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -1 /path/to/sample1_1.fastq.gz \
    -2 /path/to/sample1_2.fastq.gz \
    -p 4 \
    -o <output> \
    --resume

python3 mcclintock.py \
    -r test/sacCer2.fasta \
    -c test/sac_cer_TE_seqs.fasta \
    -1 /path/to/sample2_1.fastq.gz \
    -2 /path/to/sample2_2.fastq.gz \
    -p 4 \
    -o <output> \
    --resume

## etc ##
```
* Individual samples can be run in a serial manner as shown in the example above, or run in parallel, such as through separate jobs on a HPC cluster.

## <a name="citation"></a> Citation
To cite McClintock 1, the general TE detector meta-pipeline concept, or the single synthetic insertion simulation framework, please use: Nelson, M.G., R.S. Linheiro & C.M. Bergman (2017) McClintock: An integrated pipeline for detecting transposable element insertions in whole genome shotgun sequencing data. [G3. 7:2763-2778](https://academic.oup.com/g3journal/article/7/8/2763/6031520).

To cite McClintock 2 or the reproducible simulation system, please use: Chen, J., P.J. Basting, S. Han, D.J. Garfinkel & C.M. Bergman (2023) Reproducible evaluation of transposable element detectors with McClintock 2 guides accurate inference of Ty insertion patterns in yeast [Mobile DNA 14:8](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-023-00296-4).

## <a name="license"></a> License
------
Copyright 2014-2025 Preston Basting, Jingxuan Chen, Shunhua Han, Michael G. Nelson, and Casey M. Bergman

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
