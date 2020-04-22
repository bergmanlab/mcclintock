import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    print("installing ngs_te_mapper...")

    install_path = snakemake.config['paths']['install']+"/tools/"

    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5)

    if not download_success:
        print("temp download failed... retrying...")
        download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5)
        if not download_success:
            print("temp second download attempt failed... exiting...")
            print("try running --install with --clean for clean installation")
            sys.exit(1)

    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["mv", snakemake.config['paths']['install']+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9"):
        command = ["mv", install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9/"+f, install_path+"ngs_te_mapper"]
        mccutils.run_command(command, log=snakemake.params.log)  


    mccutils.remove(install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9")
    mccutils.remove(snakemake.params.zipfile)


    print("ngs_te_mapper installation complete")

if __name__ == "__main__":                
    main()