import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    print("installing popoolationTE...")

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

    for f in os.listdir(snakemake.config['paths']['install']+"popoolationte"):
        command = ["mv", snakemake.config['paths']['install']+"popoolationte/"+f, install_path+"popoolationte"]
        mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch", "-i", snakemake.params.patch1, install_path+"popoolationte/Modules/TEInsertUtility.pm"]
    command = ["patch", "-i", snakemake.params.patch2, install_path+"popoolationte/Modules/TEInsert.pm"]    
    command = ["patch", "-i", snakemake.params.patch3, install_path+"popoolationte/samro.pl"]
    command = ["patch", "-i", snakemake.params.patch4, install_path+"popoolationte/identify-te-insertsites.pl"]

    command = ["rm", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["rm", "-r", snakemake.config['paths']['install']+"popoolationte"]
    mccutils.run_command(command, log=snakemake.params.log)

    print("popoolationTE installation complete")

if __name__ == "__main__":                
    main()