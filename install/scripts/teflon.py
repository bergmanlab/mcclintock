import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    #set up installation path and variable names for component method
    install_path = snakemake.config['paths']['install']+"tools/"

    raw_name="TEFLoN-3e2d67886b70644fd1f7d79263b3c8dbed639e46"
    method_name = "teflon"

    #download component method source code and check for integrity
    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print(method_name+" download failed... exiting...")
        sys.exit(1)

    #unpack component method source code
    mccutils.remove(snakemake.config['paths']['install']+raw_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    #move component method source code directory into tools directory
    mccutils.remove(install_path+raw_name)
    command = ["mv", snakemake.config['paths']['install']+raw_name, install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    #move component method source code files into component method directory
    mccutils.remove(install_path+method_name)
    mccutils.mkdir(install_path+method_name)
    for f in os.listdir(install_path+raw_name):
        command = ["mv", install_path+raw_name+"/"+f, install_path+method_name]
        mccutils.run_command(command, log=snakemake.params.log)

    #patch component method source code files
    command = ["patch","-i", snakemake.params.pseudo2refConvert_patch, install_path+method_name+"/teflon_scripts/pseudo2refConvert.py"]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch","-i", snakemake.params.teflon_patch, install_path+method_name+"/teflon.v0.4.py"]
    mccutils.run_command(command, log=snakemake.params.log)

    #remove download files
    mccutils.remove(install_path+raw_name)
    mccutils.remove(snakemake.params.zipfile)

    # write version to file
    with open(install_path+method_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()