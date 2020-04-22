import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    print("<PROCESSING> making PopoolationTE reference fasta...")
    command = ["cat", snakemake.input[0], snakemake.input[1], snakemake.input[2]]
    mccutils.run_command_stdout(command, snakemake.output[0])
        

if __name__ == "__main__":                
    main()