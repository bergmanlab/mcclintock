# Reproducible and Automated McClintock Simulation Framework

McClintock provides several scripts in snakemake to simulate genomes with a single synthetic TE insertion, generate a synthetic whole genome sequence (WGS) paired-end dataset, submit the synthetic WGS dataset to the main McClintock system, and summarize component method performance.


## Install McClintock 2.0 simulation system
* Install [miniconda](https://github.com/bergmanlab/mcclintock/tree/master#-software-dependencies) and [McClintock system](https://github.com/bergmanlab/mcclintock/tree/master#-installing-mcclintock) if not done previously.

* Install the conda environments for simulation
```bash
conda activate base
cd /path/to/mcclintock
# environment for generating synthetic data and running McClintock workflow
mamba env create -f simulation/envs/mcc_sim.yml --name mcc_sim 
# environment for post-simulation data analysis
mamba env create -f simulation/envs/mcc_analysis.yml --name mcc_analysis 
```

## Quick Start
- Sample script to submit the entire evaluation framework as a SLURM cluster job:
```bash
#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=snakemake_test
#SBATCH --nodes=1
#SBATCH --tasks-per-node=10
#SBATCH --mem=40G
#SBATCH --time=100:00:00
#SBATCH --mail-user=<email.address>
#SBATCH --mail-type=BEGIN,END,FAIL

# absolute path to mcclintock repo
mcc_dir="/home/${USER}/mcclintock"
# absolute path to your expected output directory
test_dir="/scratch/${USER}/snakemake_test"
mkdir -p $test_dir
cd $test_dir

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate mcc_sim 
mkdir -p logs/ config/
# absolute path to config file for the entire simulation framework
config=${test_dir}/config/test.json

# create json config for entire simulation framework
python ${mcc_dir}/simulation/make_snakemake_config.py \
    --family "TY1" "TY2" "TY3" "TY4" \
    --bed "${mcc_dir}/simulation/yeast/200-195_tRNA_targets.bed" \
        "${mcc_dir}/simulation/yeast/200-195_tRNA_targets.bed" \
        "${mcc_dir}/simulation/yeast/17-12_tRNA_targets.bed" \
        "${mcc_dir}/simulation/yeast/200-195_tRNA_targets.bed" \
    --tsd 5 5 5 5 \
    --mcc ${mcc_dir} \
    --out ${config} \
    --outdir ${test_dir} \
    --reference "${mcc_dir}/test/sacCer2.fasta" \
    --consensus "${mcc_dir}/test/sac_cer_TE_seqs.fasta" \
    --locations "${mcc_dir}/test/reference_TE_locations.gff" \
    --taxonomy "${mcc_dir}/test/sac_cer_te_families.tsv"

# submit jobs to run simulation for each replicate
## needs to modify <email.address>
snakemake \
    --jobs 200 \
    --use-conda \
    --latency-wait 60 \
    --keep-going \
    --restart-times 3 \
    --snakefile "${mcc_dir}/simulation/Snakefile_sim" \
    --configfile ${config} \
    --cluster '
        sbatch \
            --partition=batch \
            --nodes=1 \
            --tasks-per-node={threads} \
            --mem={resources.mem} \
            --time=72:00:00 \
            --mail-user=<email.address> \
            --mail-type=FAIL \
            --output=logs/{wildcards.cov}x_{wildcards.strand}_rep{wildcards.rep}.o \
            --error=logs/{wildcards.cov}x_{wildcards.strand}_rep{wildcards.rep}.e'

# post-simulation analysis
snakemake \
    --cores $SLURM_NTASKS \
    --use-conda \
    --configfile ${config} \
    --snakefile "${mcc_dir}/simulation/Snakefile_analysis"
```

## Usage of evaluation framework
- The sample script of evaluation framework consists of three sections:
  - `simulation/make_snakemake_config.py`: Python script used to generate a configuration for simulation.
  - `simulation/Snakefile_sim`: Snakefile that creates synthetic insertions, generates synthetic WGS reads, and submits numbers of McClintock replicate jobs to the computing cluster.
  - `simulation/Snakefile_analysis`: Snakefile that runs a preliminary analysis on simulation results.

### 1: Creating configuration for simulation
```
usage: make_snakemake_config.py [-h] --family FAMILY [FAMILY ...] --bed BED
                                [BED ...] --mcc MCC --out OUT --outdir OUTDIR
                                --reference REFERENCE --consensus CONSENSUS
                                --locations LOCATIONS --taxonomy TAXONOMY
                                [--tsd TSD [TSD ...]]
                                [--covrange COVRANGE [COVRANGE ...]]
                                [--numrep NUMREP]
                                [--strand STRAND [STRAND ...]]
                                [--length LENGTH] [--insert INSERT]
                                [--error ERROR]
                                [--keep_intermediate KEEP_INTERMEDIATE]
                                [--runid RUNID] [--sim SIM] [--single]
                                [--threads THREADS] [--memory MEMORY]
                                [--exclude EXCLUDE]

Create config file in json format for snakemake pipeline that runs simulation
framework.

optional arguments:
  -h, --help            show this help message and exit

Required:
  --family FAMILY [FAMILY ...]
                        List of TE families. Required.
  --bed BED [BED ...]   List of candidate insertion region files in BED format
                        for each TE family, respectively. Or one bed file for
                        all TE families. Required.
  --mcc MCC             Path to local McClintock repo. Required.
  --out OUT             File name of the output config file in json. Required.
  --outdir OUTDIR       Output directory for all simulation results. Required.
  --reference REFERENCE
                        [McClintock option] A reference genome sequence in
                        fasta format. Required.
  --consensus CONSENSUS
                        [McClintock option] The consensus sequences of the TEs
                        for the species. Required.
  --locations LOCATIONS
                        [McClintock option] The locations of known TEs in the
                        reference genome in GFF 3 format. Required.
  --taxonomy TAXONOMY   [McClintock option] A tab delimited file with one
                        entry per ID. Required.

Optional:
  --tsd TSD [TSD ...]   Integers for TSD length of each TE family in bp.
                        [default = 5]
  --covrange COVRANGE [COVRANGE ...]
                        List of ingeters for simulated WGS coverage range.
                        [default = 3 6 12 25 50 100]
  --numrep NUMREP       Integer of replicate counts for each coverage and each
                        strand. [default = 30]
  --strand STRAND [STRAND ...]
                        List of strands for simulation. 'forward' and
                        'reverse' are allowed. [default = 'forward' 'reverse']

Optional simulation parameters:
  --length LENGTH       The read length of the simulated reads [default = 101]
  --insert INSERT       The median insert size of the simulated reads [default
                        = 300]
  --error ERROR         The base error rate for the simulated reads [default =
                        0.01]
  --keep_intermediate KEEP_INTERMEDIATE
                        This option determines which intermediate files are
                        preserved after McClintock completes (options:
                        minimal, general, methods, <list,of,methods>,
                        all)[default = 'minimal']
  --runid RUNID         (not recommended) a string to prepend to output files
                        so that multiple runs can be run at the same time
                        without file name clashes
  --sim SIM             Short read simulator to use (options: wgsim,art)
                        [default = 'art']
  --single              Runs the simulation in single end mode.

Optional resources:
  --threads THREADS     The number of processors to use for individual cluster
                        jobs. [default = 4]
  --memory MEMORY       The number of memory in 'G' to use for individual
                        cluster jobs. [default = '20G']

Optional analysis parameters:
  --exclude EXCLUDE     BED file of regions in which predictions will be
                        excluded from counts (ex. low recombination regions),
                        for analysis script.
```
#### `simulation/make_snakemake_config.py` inputs
- Required:
  - `--reference`: A reference genome sequence in FASTA format. 
    - This is the same required input file to run McClintock.
    - A sample for yeast could be found at `test/sacCer2.fasta`.
  - `--consensus`: The consensus sequences of the TEs for the species in FASTA format.
    - This is the same required input file to run McClintock.
    - A sample for yeast could be found at `test/sac_cer_TE_seqs.fasta`.
  - `--locations`: The locations of known TEs in the reference genome in GFF3 format.
    - This is the same optional input file to run McClintock but required for simulation framework.
    - Locations GFF3 could be generated by running `python3 <mcc_dir>/mcclintock.py --make_annotations -r <reference> -c <consensus> -o <output_folder>` 
      - `<output_folder>/<prefix_of_ref>/reference_te_locations/inrefTEs.gff`
    - A sample for yeast could be found at `test/reference_TE_locations.gff`
  - `--taxonomy`: A tab delimited file with one entry per ID.
    - This is the same optional input file to run McClintock but required for simulation framework.
    - Taxonomy TSV could be generated by running `python3 <mcc_dir>/mcclintock.py --make_annotations -r <reference> -c <consensus> -o <output_folder>`
      - `<output_folder>/<prefix_of_ref>/te_taxonomy/taxonomy.tsv`
    - A sample for yeast could be found at `test/sac_cer_te_families.tsv`
  - `--family`: Space delimited list of TE families in the species. TE names must be consistent with the consensus sequences (`--consensus`).
  - `--bed`: Space delimited list of candidate insertion region files in BED format for each TE family, respectively. Or one bed file for all TE families.
  - `--mcc`: Absolute path to local McClintock repo.
  - `--out`: File name for the output config file in JSON.
  - `--outdir`: Output directory for the entire simulation framework results.
- Optional:
  - General optional parameters:
    - `--tsd`: Space delimited list of TSD length for each TE family in bp. Or one integer for all TE families. 
      - Default: 5 bp for all TE families.
    - `--covrange`: Space delimited list of simulated WGS data coverage range.
      - Default: 3 6 12 25 50 100.
    - `--numrep`: Number of replicates for each simulated coverage and each strand.
      - Default: 30.
    - `--strand`: Strands to create synthetic insertion. 'forward' and 'reverse' are allowed.
      - Default: 'forward' 'reverse'.
  - Parameters to run simulation:
    - `--length`: The read length of the simulated WGS reads in bp.
      - Default: 101.
    - `--insert`: The median insert size of the simulated WGS reads in bp.
      - Default: 300.
    - `--error`: The base error rate for the simulated WGS reads. Only useful when using wgsim as the read simulator.
      - Default: 0.01.
    - `--keep_intermediate`: Intermediate files to preserve after McClintock completes. This is helpful when users would keep large intermediate files (eg. SAM files) for debug.
      - This is the same optional input file to run McClintock.
      - Options: `minimal, general, methods, <list,of,methods>, all`.
      - Default: minimal.
    - `--sim`: Specify the short-read simulator to use.
      - Options: wgsim,art.
      - Default: art.
    - `--single`: If specified, runs read simulator and McClintock in single end mode.
      - Default: Not specified, run in paired-end mode.
  - Parameters for computing resources:
    - `--threads`: Number of processors to run simulater and McClintock for each replicate.
      - Default: 4.
    - `--memory`: The amount of memory to run simulater and McClintock for each replicate in G.
      - Default: 20G.
  - Parameters for analysis:
    - `--exclude`: Regions in which predictions will be excluded from counts (ex. low recombination regions) in BED format.
- A sample config file in json could be found at `simulation/cfg_sim.json`.

### 2: Submitting individual cluster jobs for simulation replicates
- In this section, we implemented a "single synthetic insertion" simulation framework to evaluate the performance of McClintock component methods.
  - Firstly, a synthetic genome comprised of a single non-reference TE insertion and its corresponding target site duplication (TSD) is created.
  - Then, synthetic paired-end WGS dataset are generated with read simulator.
  - Finally, McClintock meta-pipeline is executed with simulated WGS data as input.

![singleinssim](https://user-images.githubusercontent.com/44645406/204376432-a1d617d4-64b7-4a3f-aa06-5fe95a86e89f.jpeg)


- Each replicate for simulation framework would be submitted as a separate cluster job that could run in parallel.
- Run Snakefile `simulation/Snakefile_sim` along with the config file:
```
conda activate mcc_sim 
snakemake \
    --jobs 200 \
    --use-conda \
    --latency-wait 60 \
    --keep-going \
    --restart-times 3 \
    --snakefile "${mcc_dir}/simulation/Snakefile_sim" \
    --configfile ${config} \
    --cluster '
        sbatch \
            --partition=batch \
            --nodes=1 \
            --tasks-per-node={threads} \
            --mem={resources.mem} \
            --time=72:00:00 \
            --mail-user=<email.address> \
            --mail-type=FAIL \
            --output=logs/{wildcards.cov}x_{wildcards.strand}_rep{wildcards.rep}.o \
            --error=logs/{wildcards.cov}x_{wildcards.strand}_rep{wildcards.rep}.e'
```
- Parameters to run `snakemake`
  - `--jobs`: The number of jobs to be submitted at the same time. Please modify this parameter according to available amount of computing resources on your cluster.
  - `--use-conda`: Required. Use conda environment in snakemake rule.
  - `--latency-wait`: Wait given seconds if an output file of a job is not present after the job finished. Please modify this parameter according to latency of your filesystem.
  - `--keep-going`: Required. Keep running even though there are failed jobs.
  - `--restart-times`: Number of times to re-run if job failed. Recommended.
  - `--snakefile`: Required. Specify the absolute path to `simulation/Snakefile_sim`.
  - `--configfile`: Required. Specify the absolute path to the configuration file.
  - `--cluster`: Required. Provide headers for job submission. *Note that if `--mail-user` and `--mail-type` are provided, you may get tons of emails for each replicate.*
  - Other arguments to run snakemake may be applicable. Please check [snakemake tutorial](https://snakemake.readthedocs.io/en/v7.16.0/index.html).
- External scripts called by `mcclintock/simulation/Snakefile_sim`:
  - `mcclintock/simulation/make_sim_config.py`
  - `mcclintock/simulation/mcclintock_simulation_snk.py`

### 3: Preliminary anslysis on simulation results
- Run Snakefile `simulation/Snakefile_analysis` along with the config file:
```
conda activate mcc_sim 
snakemake \
    --cores $SLURM_NTASKS \
    --use-conda \
    --configfile ${config} \
    --snakefile "${mcc_dir}/simulation/Snakefile_analysis"
```
- Output files:
  - `<outdir>/<cov>/summary/combined_metrics.csv`: Summary matrix for precision and recall, as well as raw numbers of true-positive (TP), false-positive (FP) and false-negative (FN).
  - `<outdir>/<cov>/summary/forward_metrics.csv` and `<outdir>/<cov>/summary/reverse_metrics.csv`: Summary matrix for precision and recall for each strand separately.
  - `<outdir>/<cov>/summary/forward_summary.csv` and `<outdir>/<cov>/summary/reverse_summary.csv`: Average numbers of reference and non-reference TE predictions across replicates. Like table 4 shown in [Nelson, Linheiro and Bergman (2017) G3 7:2763-2778](https://academic.oup.com/g3journal/article/7/8/2763/6031520).
  - `<outdir>/r_vis`: Visualization of simulation results
    - `<outdir>/r_vis/sim_<Precision/Recall>_mainfig.pdf` and `<outdir>/r_vis/sim_<Precision/Recall>_allcovs.pdf`: Precision and recall curves for 2 window sizes and all windows, respectively.
    - `<outdir>/r_vis/simdensity_<cov>x_<TE_name>.pdf`: Positional accuracy. Distribution of non-reference predictions relative to true location of component methods across fold-coverage and TE families.
    - `<outdir>/r_vis/simtsd_<cov>x_<TE_name>.pdf`: TSD length disbutions for component methods using split-read evidence.
- External scripts called by `mcclintock/simulation/Snakefile_analysis`:
  - `mcclintock/simulation/mcclintock_simulation_analysis.py`
  - `mcclintock/simulation/sim_r_vis_precision_recall.R`
  - `mcclintock/simulation/sim_r_vis_for_cov.R`
