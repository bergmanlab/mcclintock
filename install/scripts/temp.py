import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils
import subprocess

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("temp download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)

    extracted_file_name = "TEMP-4f67e1da836721a9f0999efa52e1e648fedb75fc"

    mccutils.remove(snakemake.config['paths']['install']+extracted_file_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+extracted_file_name)
    command = ["mv", snakemake.config['paths']['install']+extracted_file_name, install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"temp")
    mccutils.mkdir(install_path+"temp")
    for f in os.listdir(install_path+extracted_file_name):
        command = ["mv", install_path+extracted_file_name+"/"+f, install_path+"temp"]
        mccutils.run_command(command, log=snakemake.params.log)


    command = ["patch", "-i", snakemake.params.patch, install_path+"temp/scripts/TEMP_Absence.sh"]
    mccutils.run_command(command,log=snakemake.params.log) 

    mccutils.remove(install_path+extracted_file_name)
    mccutils.remove(snakemake.params.zipfile)

    # write version to file
    with open(snakemake.config['paths']['install']+"/tools/temp/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()