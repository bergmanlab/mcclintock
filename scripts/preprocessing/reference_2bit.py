import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    log = snakemake.params.log
    print("<PROCESSING> creating 2bit file from reference genome fasta...log:"+log)
    command = ["faToTwoBit", snakemake.input[0], snakemake.output[0]]
    mccutils.run_command(command, log=log)
        

if __name__ == "__main__":                
    main()