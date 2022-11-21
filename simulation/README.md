# McClintock 2.0 Simulation System

McClintock 2.0 provides several scripts to simulate genomes with a single synthetic TE insertion, generate a synthetic whole genome sequence (WGS) paired-end dataset, submit the synthetic WGS dataset to the main McClintock system, and summarize component method performance.


## Install the McClintock 2.0 simulation system
* Install miniconda (if not done previously)
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda # silent mode
echo "export PATH=\$PATH:\$HOME/miniconda/bin" >> $HOME/.bashrc # add to .bashrc
source $HOME/.bashrc
conda init
# close and re-open terminal
conda update conda
conda install -c conda-forge mamba
```

* Clone and install McClintock 2.0 (if not done previously)
```bash
cd ~
git clone git@github.com:bergmanlab/mcclintock.git
cd mcclintock

# install and activate the main McClintock environment
mamba env create -f install/envs/mcclintock.yml --name mcclintock
CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate mcclintock

# install McClintock component methods and environments
python3 mcclintock.py --install
```

* Install the conda environments for simulation
```bash
# environment for generating synthetic data and runing McClintock workflow
mamba env create -f /path/to/mcclintock/simulation/envs/mcc_sim.yml --name mcc_sim
# environment for post-simulation data analysis
mamba env create -f /path/to/mcclintock/simulation/envs/mcc_analysis.yml --name mcc_analysis
```

## Quick Start
- Submit the following script as a cluster job:

```bash
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
            --mail-type=END,FAIL \
            --output=logs/{wildcards.cov}x_{wildcards.strand}_rep{wildcards.rep}.o \
            --error=logs/{wildcards.cov}x_{wildcards.strand}_rep{wildcards.rep}.e'

# post-simulation analysis
snakemake \
    --cores $SLURM_NTASKS \
    --use-conda \
    --configfile ${config} \
    --snakefile "${mcc_dir}/simulation/Snakefile_analysis"
```
