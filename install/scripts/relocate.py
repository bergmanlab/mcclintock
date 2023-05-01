import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils
import subprocess

def main():
    #set up installation path and variable names for component method
    install_path = snakemake.config['paths']['install']+"/tools/"

    raw_name = "RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b"
    method_name = "relocate"

    #download component method source code and check for integrity
    #note: snakemake.params are found in /install/Snakefile
    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print(method_name+" download failed... exiting...")
        sys.exit(1)

    #unpack component method source code
    mccutils.remove(snakemake.config['paths']['install']+raw_name)
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    #move component method source code directory into tools directory
    mccutils.remove(install_path+raw_name)
    command = ["mv", snakemake.config['paths']['install']+raw_name, install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    #move component method source code files into component method directory
    mccutils.remove(install_path+method_name)
    mccutils.mkdir(install_path+method_name)
    for f in os.listdir(install_path+raw_name):
        command = ["mv", install_path+raw_name+"/"+f, install_path+method_name]
        mccutils.run_command(command, log=snakemake.params.log)

    #patch component method source code files
    command = ["patch", "-i", snakemake.params.patch, install_path+method_name+"/scripts/relocaTE_insertionFinder.pl"]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+raw_name)
    mccutils.remove(snakemake.params.zipfile)

    output = subprocess.Popen(["which", "perl"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    perl_path = output.stdout.read()
    perl_path = perl_path.decode()

    #update version of perl used in source code files
    for f in os.listdir(install_path+method_name+"/scripts/"):
        if "pl" == f.split(".")[-1]:
            with open(install_path+"tmp","w") as tmp:
                with open(install_path+method_name+"/scripts/"+f, "r") as script:
                    for line in script:
                        if "#!/usr/bin/perl" in line:
                            # line = "#!"+perl_path
                            line = "#!/usr/bin/env perl\n"
                        elif "defined @" in line:
                            line = line.replace("defined @", "@")
                        
                        elif "$scripts/" in line and "perl" not in line and "relocaTE.pl" in f:
                            line = line.replace("$scripts/", "perl $scripts/")
                        
                        tmp.write(line)
            
            mccutils.run_command(["mv", install_path+"tmp", install_path+method_name+"/scripts/"+f])

    #write version to file
    with open(install_path+method_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()