import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():

    command = ["wget", "--no-check-certificate", "https://downloads.sourceforge.net/project/te-locate/TE-locate.tar", "-O", snakemake.params.tar]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["tar", "-xvf", snakemake.params.tar, "-C", snakemake.config['paths']['install']+"/tools/te-locate/"]
    mccutils.run_command(command, log=snakemake.params.log)

if __name__ == "__main__":                
    main()