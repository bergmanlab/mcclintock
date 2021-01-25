import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    raw_name = "jitterbug-b6b3f9c7ee4af042d4410137269f174a1399b752"
    mod_name = "jitterbug"

    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print(mod_name+" download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)

    mccutils.remove(snakemake.config['paths']['install']+raw_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+raw_name)
    command = ["mv", snakemake.config['paths']['install']+raw_name, install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+mod_name)
    mccutils.mkdir(install_path+mod_name)

    for f in os.listdir(install_path+raw_name):
        command = ["mv", install_path+raw_name+"/"+f, install_path+mod_name]
        mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch", "-i", snakemake.params.patch, install_path+mod_name+"/jitterbug.py"]
    mccutils.run_command(command, log=snakemake.params.log)
    mccutils.remove(install_path+raw_name)
    mccutils.remove(snakemake.params.zipfile)


    # write version to file
    with open(install_path+mod_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()