import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    print("installing popoolationTE2...")

    download_success = mccutils.download(snakemake.params.url, snakemake.output[0], md5=snakemake.params.md5)
    if not download_success:
        print("popoolationTE2 download failed... retrying...")
        download_success = mccutils.download(snakemake.params.url, snakemake.output[0], md5=snakemake.params.md5)
        if not download_success:
            print("second popoolationTE2 download attempt failed... exiting...")
            print("try running --install with --clean for clean installation")
            sys.exit(1)

    print("popoolationTE2 installation complete")

if __name__ == "__main__":                
    main()