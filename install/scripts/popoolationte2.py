import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    download_success = mccutils.download(snakemake.params.url, snakemake.output[0], md5=snakemake.params.md5, max_attempts=3)
    if not download_success:
        print("popoolationTE2 download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)

if __name__ == "__main__":                
    main()