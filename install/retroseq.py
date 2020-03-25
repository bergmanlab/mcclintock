import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    mccutils.mkdir(install_path)

    command = ["wget", "--no-check-certificate", "https://github.com/tk2/RetroSeq/archive/700d4f76a3b996686652866f2b81fefc6f0241e0.zip", "-O", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["mv", snakemake.config['paths']['install']+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0"):
        command = ["mv", install_path+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0/"+f, install_path+"retroseq"]
        mccutils.run_command(command, log=snakemake.params.log)
    
    command = ["rm", "-r", install_path+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0"]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["rm", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)


if __name__ == "__main__":                
    main()