import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"

    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("ngs_te_mapper download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)

    mccutils.remove(snakemake.config['paths']['install']+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9")
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9")
    command = ["mv", snakemake.config['paths']['install']+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"ngs_te_mapper")
    mccutils.mkdir(install_path+"ngs_te_mapper")
    for f in os.listdir(install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9"):
        command = ["mv", install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9/"+f, install_path+"ngs_te_mapper"]
        mccutils.run_command(command, log=snakemake.params.log)  


    mccutils.remove(install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9")
    mccutils.remove(snakemake.params.zipfile)

if __name__ == "__main__":                
    main()