import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    print("<PROCESSING> sorting SAM file for compatibility with TE-locate...")
    command = ["sort", "-S", snakemake.config['args']['mem']+"G", "--temporary-directory="+snakemake.config['args']['out']+"/tmp", snakemake.input[0]]
    mccutils.run_command_stdout(command, snakemake.output[0])
        

if __name__ == "__main__":                
    main()