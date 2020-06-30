import sys
import os
import subprocess
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("RelocaTE2 download failed... exiting...")
        print("try running --install with --clean for clean installation")
        sys.exit(1)

    mccutils.remove(snakemake.config['paths']['install']+"RelocaTE2-2.0.1")
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"RelocaTE2-2.0.1")
    command = ["mv", snakemake.config['paths']['install']+"RelocaTE2-2.0.1", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"relocate2")
    mccutils.mkdir(install_path+"relocate2")
    for f in os.listdir(install_path+"RelocaTE2-2.0.1"):
        command = ["mv", install_path+"RelocaTE2-2.0.1/"+f, install_path+"relocate2"]
        mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"RelocaTE2-2.0.1")
    mccutils.remove(snakemake.params.zipfile)

    tools = [ "bwa", "bowtie2", "bowtie2_build", "blat", "samtools", "bedtools", "seqtk" ]

    with open(install_path+"relocate2/CONFIG","w") as config_file:
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


if __name__ == "__main__":                
    main()