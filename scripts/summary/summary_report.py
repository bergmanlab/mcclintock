import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    out_files = snakemake.input.out_files
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2
    ref = snakemake.input.ref
    bam = snakemake.input.bam
    flagstat = snakemake.input.flagstat
    median_insert_size = snakemake.input.median_insert_size

    
    
    mccutils.run_command(["touch", snakemake.output[0]])

if __name__ == "__main__":                
    main()