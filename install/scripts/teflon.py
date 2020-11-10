import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"

    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("teflon download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)
    
    mccutils.remove(snakemake.config['paths']['install']+"TEFLoN-9eca0152f3dd9dc6c44787a30d590f7e321b442c")
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"TEFLoN-9eca0152f3dd9dc6c44787a30d590f7e321b442c")
    command = ["mv", snakemake.config['paths']['install']+"TEFLoN-9eca0152f3dd9dc6c44787a30d590f7e321b442c", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"teflon")
    mccutils.mkdir(install_path+"teflon")

    for f in os.listdir(install_path+"TEFLoN-9eca0152f3dd9dc6c44787a30d590f7e321b442c"):
        command = ["mv", install_path+"TEFLoN-9eca0152f3dd9dc6c44787a30d590f7e321b442c/"+f, install_path+"teflon"]
        mccutils.run_command(command, log=snakemake.params.log)


    command = ["patch","-i", snakemake.params.pseudo2refConvert_patch, install_path+"teflon/teflon_scripts/pseudo2refConvert.py"]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["patch","-i", snakemake.params.teflon_patch, install_path+"teflon/teflon.v0.4.py"]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"TEFLoN-9eca0152f3dd9dc6c44787a30d590f7e321b442c")
    mccutils.remove(snakemake.params.zipfile)

if __name__ == "__main__":                
    main()