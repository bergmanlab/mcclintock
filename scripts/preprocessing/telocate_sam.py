import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    log = snakemake.params.log
    mccutils.log("processing","sorting SAM file for compatibility with TE-locate", log=log)
    command = ["sort", "-S", "1G", "--temporary-directory="+snakemake.config['args']['out']+"/tmp", snakemake.input[0]]
    mccutils.run_command_stdout(command, snakemake.output[0], log=log)
    mccutils.log("processing","TE-locate SAM created")
        

if __name__ == "__main__":                
    main()