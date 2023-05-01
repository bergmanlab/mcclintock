import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils
import subprocess

def main():
    #set up installation path and variable names for component method
    install_path = snakemake.config['paths']['install']+"tools/"

    raw_name = "RelocaTE2-2.0.1"
    method_name = "relocate2"

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

    mccutils.remove(install_path+raw_name)
    mccutils.remove(snakemake.params.zipfile)

    tools = [ "bwa", "bowtie2", "bowtie2_build", "blat", "samtools", "bedtools", "seqtk" ]

    #update paths to executables in relocate2 config file
    with open(install_path+method_name"/CONFIG","w") as config_file:
        config_file.write("#tools\n")
        for tool in tools:
            if tool == "bowtie2_build":
                tool_bin = "bowtie2-build"
            else:
                tool_bin = tool
            output = subprocess.Popen(["which", tool_bin], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            path = output.stdout.read()
            path = path.decode()
            line = tool+"="+path
            config_file.write(line)

    #write version to file
    with open(install_path+method_name+"/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()