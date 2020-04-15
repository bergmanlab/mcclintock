import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils

def main():
    relocate_gff = snakemake.input.relocate_gff
    te_gff = snakemake.input.te_gff

    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name

    mccutils.run_command(["touch", snakemake.output[0]])


if __name__ == "__main__":                
    main()