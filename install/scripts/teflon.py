import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"

    raw_name="TEFLoN-3e2d67886b70644fd1f7d79263b3c8dbed639e46"
    method_name = "teflon"

    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print(method_name+" download failed... exiting...")
        sys.exit(1)
    
    mccutils.remove(snakemake.config['paths']['install']+raw_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+raw_name)
    command = ["mv", snakemake.config['paths']['install']+raw_name, install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"teflon")
    mccutils.mkdir(install_path+"teflon")

    for f in os.listdir(install_path+raw_name):
        command = ["mv", install_path+raw_name+"/"+f, install_path+"teflon"]
        mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch","-i", snakemake.params.pseudo2refConvert_patch, install_path+method_name+"/teflon_scripts/pseudo2refConvert.py"]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch","-i", snakemake.params.teflon_patch, install_path+method_name+"/teflon.v0.4.py"]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+raw_name)
    mccutils.remove(snakemake.params.zipfile)

    # write version to file
    with open(install_path+"/"+method_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()