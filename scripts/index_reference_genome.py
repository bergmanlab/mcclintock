import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    log = snakemake.params.log
    print("<PROCESSING> making samtools and bwa index files for reference fasta")
    mccutils.run_command(["samtools", "faidx", snakemake.input.ref],log=log)
    mccutils.run_command(["bwa", "index", snakemake.input.ref], log=log)


if __name__ == "__main__":                
    main()