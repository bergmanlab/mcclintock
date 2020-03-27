import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    command = ["faToTwoBit", snakemake.input[0], snakemake.output[0]]
    mccutils.run_command(command)
        

if __name__ == "__main__":                
    main()