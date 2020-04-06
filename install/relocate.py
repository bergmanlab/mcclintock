import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    print("installing relocaTE...")

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

    command = ["mv", snakemake.config['paths']['install']+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b"):
        command = ["mv", install_path+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b/"+f, install_path+"relocate"]
        mccutils.run_command(command, log=snakemake.params.log)


    command = ["patch", "-i", snakemake.params.patch, install_path+"relocate/scripts/relocaTE_insertionFinder.pl"]    

    mccutils.remove(install_path+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b")
    mccutils.remove(snakemake.params.zipfile)

    print("relocaTE installation complete")

if __name__ == "__main__":                
    main()