import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    print("installing RetroSeq...")

    install_path = snakemake.config['paths']['install']+"/tools/"

    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5)

    if not download_success:
        print("retroseq download failed... retrying...")
        download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5)
        if not download_success:
            print("retroseq second download attempt failed... exiting...")
            print("try running --install with --clean for clean installation")
            sys.exit(1)

    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["mv", snakemake.config['paths']['install']+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0"):
        command = ["mv", install_path+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0/"+f, install_path+"retroseq"]
        mccutils.run_command(command, log=snakemake.params.log)
    
    mccutils.remove(install_path+"RetroSeq-700d4f76a3b996686652866f2b81fefc6f0241e0")
    mccutils.remove(snakemake.params.zipfile)

    print("RetroSeq installation complete")

if __name__ == "__main__":                
    main()