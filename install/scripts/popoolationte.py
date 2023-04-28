import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"

    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("popoolationte download failed... exiting...")
        sys.exit(1)

    mccutils.remove(snakemake.config['paths']['install']+"popoolationte")
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"popoolationte")
    mccutils.mkdir(install_path+"popoolationte")
    for f in os.listdir(snakemake.config['paths']['install']+"popoolationte"):
        command = ["mv", snakemake.config['paths']['install']+"popoolationte/"+f, install_path+"popoolationte"]
        mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch", "-i", snakemake.params.patch1, install_path+"popoolationte/Modules/TEInsertUtility.pm"]
    mccutils.run_command(command, log=snakemake.params.log)
    command = ["patch", "-i", snakemake.params.patch2, install_path+"popoolationte/Modules/TEInsert.pm"]
    mccutils.run_command(command, log=snakemake.params.log)
    command = ["patch", "-i", snakemake.params.patch3, install_path+"popoolationte/samro.pl"]
    mccutils.run_command(command, log=snakemake.params.log)
    command = ["patch", "-i", snakemake.params.patch4, install_path+"popoolationte/identify-te-insertsites.pl"]
    mccutils.run_command(command, log=snakemake.params.log)


    mccutils.remove(snakemake.params.zipfile)
    mccutils.remove(snakemake.config['paths']['install']+"popoolationte")

    # write version to file
    with open(snakemake.config['paths']['install']+"/tools/popoolationte/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()