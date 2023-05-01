import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    #set up installation path and variable names for component method
    install_path = snakemake.config['paths']['install']+"tools/"

    method_name = "te-locate"

    #download component method source code and check for integrity
    #note: snakemake.params are found in /install/Snakefile
    mccutils.remove(snakemake.params.tar)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.tar, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print(method_name+" download failed... exiting...")
        sys.exit(1)

    #unpack component method source code into component method directory
    command = ["tar", "-xvf", snakemake.params.tar, "-C", install_path+method_name]
    mccutils.run_command(command, log=snakemake.params.log)

    #write version to file
    with open(install_path+method_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()