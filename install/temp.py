import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import modules.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    mccutils.mkdir(install_path)

    command = ["wget", "--no-check-certificate", "https://github.com/JialiUMassWengLab/TEMP/archive/d2500b904e2020d6a1075347b398525ede5feae1.zip", "-O", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["mv", snakemake.config['paths']['install']+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1"):
        command = ["mv", install_path+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1/"+f, install_path+"temp"]
        mccutils.run_command(command, log=snakemake.params.log)


    command = ["patch", "-i", snakemake.params.patch, install_path+"temp/scripts/TEMP_Absence.sh"]    

    command = ["rm", "-r", install_path+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1"]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["rm", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

if __name__ == "__main__":                
    main()