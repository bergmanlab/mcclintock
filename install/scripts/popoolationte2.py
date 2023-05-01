import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    #set up installation path and variable names for component method
    install_path = snakemake.config['paths']['install']+"tools/"

    method_name = "popoolationte2"

    #download component method source code and check for integrity
    #note: snakemake.params are found in /install/Snakefile
    #download is .jar file, nothing to unpack
    download_success = mccutils.download(snakemake.params.url, snakemake.output[0], md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print(method_name+" download failed... exiting...")
        sys.exit(1)

    #write version to file
    with open(install_path+method_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()