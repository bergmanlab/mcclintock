import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"


    method_name = "popoolationte"

    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print(method_name+" download failed... exiting...")
        sys.exit(1)

    mccutils.remove(snakemake.config['paths']['install']+method_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+method_name)
    mccutils.mkdir(install_path+method_name)
    for f in os.listdir(snakemake.config['paths']['install']+method_name):
        command = ["mv", snakemake.config['paths']['install']+"popoolationte/"+f, install_path+method_name]
        mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch", "-i", snakemake.params.patch1, install_path+"/"+method_name+"/Modules/TEInsertUtility.pm"]
    mccutils.run_command(command, log=snakemake.params.log)
    command = ["patch", "-i", snakemake.params.patch2, install_path+"/"+method_name+"/Modules/TEInsert.pm"]
    mccutils.run_command(command, log=snakemake.params.log)
    command = ["patch", "-i", snakemake.params.patch3, install_path+"/"+method_name+"/samro.pl"]
    mccutils.run_command(command, log=snakemake.params.log)
    command = ["patch", "-i", snakemake.params.patch4, install_path+"/"+method_name+"/identify-te-insertsites.pl"]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(snakemake.params.zipfile)
    mccutils.remove(snakemake.config['paths']['install']+method_name)

    # write version to file
    with open(install_path+"/"+method_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()