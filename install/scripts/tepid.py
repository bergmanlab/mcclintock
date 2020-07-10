import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils
import subprocess

def main():
    install_path = snakemake.config['paths']['install']+"/tools/tepid/"

    # patch tepid-map
    path_file = install_path+"path.txt"
    mccutils.run_command_stdout(["which","tepid-map"],path_file)
    path = ""
    with open(path_file,"r") as inf:
        for line in inf:
            path += line
            
    path = path.replace("\n","")
    command = ["patch", "-i", snakemake.params.map_patch, path, "-o", install_path+"tepid-map"]
    mccutils.run_command(command,log=snakemake.params.log)


    # patch tepid-map-se
    path_file = install_path+"path.txt"
    mccutils.run_command_stdout(["which","tepid-map-se"],path_file)
    path = ""
    with open(path_file,"r") as inf:
        for line in inf:
            path += line

    path = path.replace("\n","")
    command = ["patch", "-i", snakemake.params.map_se_patch, path, "-o", install_path+"tepid-map-se"]
    mccutils.run_command(command,log=snakemake.params.log)

    mccutils.remove(path_file)

if __name__ == "__main__":                
    main()