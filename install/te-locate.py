import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    print("installing TE-locate...")

    download_success = mccutils.download(snakemake.params.url, snakemake.params.tar, md5=snakemake.params.md5)

    if not download_success:
        print("te-locate download failed... retrying...")
        download_success = mccutils.download(snakemake.params.url, snakemake.params.tar, md5=snakemake.params.md5)
        if not download_success:
            print("retroseq second download attempt failed... exiting...")
            print("try running --install with --clean for clean installation")
            sys.exit(1)


    command = ["tar", "-xvf", snakemake.params.tar, "-C", snakemake.config['paths']['install']+"/tools/te-locate/"]
    mccutils.run_command(command, log=snakemake.params.log)

    print("TE-locate installation complete")

if __name__ == "__main__":                
    main()