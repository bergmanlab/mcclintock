import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    mccutils.mkdir(install_path)

    command = ["wget", "--no-check-certificate", "https://github.com/bergmanlab/ngs_te_mapper/archive/fb23590200666fe66f1c417c5d5934385cb77ab9.zip", "-O", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["mv", snakemake.config['paths']['install']+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9"):
        command = ["mv", install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9/"+f, install_path+"ngs_te_mapper"]
        mccutils.run_command(command, log=snakemake.params.log)  

    command = ["rm", "-r", install_path+"ngs_te_mapper-fb23590200666fe66f1c417c5d5934385cb77ab9"]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["rm", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

if __name__ == "__main__":                
    main()