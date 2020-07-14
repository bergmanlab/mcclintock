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
        print("TEPID download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)
    
    extracted_file_name = "TEPID-ad46d65b5c41bf8a9171215d49b3ffaecdceaab0"

    mccutils.remove(snakemake.config['paths']['install']+extracted_file_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+extracted_file_name)
    command = ["mv", snakemake.config['paths']['install']+extracted_file_name, install_path+extracted_file_name]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"tepid")
    mccutils.mkdir(install_path+"tepid")
    for f in os.listdir(install_path+extracted_file_name):
        command = ["mv", install_path+extracted_file_name+"/"+f, install_path+"tepid"]
        mccutils.run_command(command, log=snakemake.params.log)  

    for f in os.listdir(install_path+"tepid/Scripts/"):
        command = ["cp", install_path+"tepid/Scripts/"+f, install_path+"tepid/"]
        mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch", "-i", snakemake.params.discover_patch, install_path+"tepid/tepid-discover"]
    mccutils.run_command(command,log=snakemake.params.log) 

    mccutils.remove(install_path+extracted_file_name)
    mccutils.remove(snakemake.params.zipfile)

if __name__ == "__main__":                
    main()