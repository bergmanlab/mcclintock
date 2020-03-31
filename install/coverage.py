import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    print("installing coverage module...")

    install_path = snakemake.config['paths']['install']+"/tools/"

    command = ["wget", "--no-check-certificate", "https://github.com/bergmanlab/samplot/archive/1de65afd22e88c5cb5122ae638e8ba4cf6f75144.zip", "-O", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["mv", snakemake.config['paths']['install']+"samplot-1de65afd22e88c5cb5122ae638e8ba4cf6f75144", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"samplot-1de65afd22e88c5cb5122ae638e8ba4cf6f75144"):
        command = ["mv", install_path+"samplot-1de65afd22e88c5cb5122ae638e8ba4cf6f75144/"+f, install_path+"coverage"]
        mccutils.run_command(command, log=snakemake.params.log)

    command = ["rm", "-r", install_path+"samplot-1de65afd22e88c5cb5122ae638e8ba4cf6f75144"]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["rm", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    print("coverage module installation complete")

if __name__ == "__main__":                
    main()