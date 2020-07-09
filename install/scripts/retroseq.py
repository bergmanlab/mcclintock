import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("retroseq download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)

    extracted_file_name = "RetroSeq-9d4f3b5270af2383f40e6e7ea1204ea718365db2"

    mccutils.remove(snakemake.config['paths']['install']+extracted_file_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+extracted_file_name)
    command = ["mv", snakemake.config['paths']['install']+extracted_file_name, install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"retroseq")
    mccutils.mkdir(install_path+"retroseq")
    for f in os.listdir(install_path+extracted_file_name):
        command = ["mv", install_path+extracted_file_name+"/"+f, install_path+"retroseq"]
        mccutils.run_command(command, log=snakemake.params.log)
    
    mccutils.remove(install_path+extracted_file_name)
    mccutils.remove(snakemake.params.zipfile)

if __name__ == "__main__":                
    main()