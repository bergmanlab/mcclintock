import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"

    raw_name = "ngs_te_mapper-f9f48996ac346ac86d57edbd00534aa1227b753e"

    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("ngs_te_mapper download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)

    mccutils.remove(snakemake.config['paths']['install']+raw_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+raw_name)
    command = ["mv", snakemake.config['paths']['install']+raw_name, install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"ngs_te_mapper")
    mccutils.mkdir(install_path+"ngs_te_mapper")
    for f in os.listdir(install_path+raw_name):
        command = ["mv", install_path+raw_name+"/"+f, install_path+"ngs_te_mapper"]
        mccutils.run_command(command, log=snakemake.params.log)  


    mccutils.remove(install_path+raw_name)
    mccutils.remove(snakemake.params.zipfile)

    # write version to file
    with open(install_path+"ngs_te_mapper/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()