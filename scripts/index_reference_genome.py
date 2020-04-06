import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    print("<PROCESSING> making samtools and bwa index files for reference fasta")
    mccutils.run_command(["samtools", "faidx", snakemake.input.ref])
    mccutils.run_command(["bwa", "index", snakemake.input.ref])


if __name__ == "__main__":                
    main()