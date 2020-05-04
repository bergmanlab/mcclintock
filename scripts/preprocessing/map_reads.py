import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    print("<PROCESSING> mapping reads to reference...log:"+snakemake.log[0])
    command = ["bwa","mem"]
    if eval(snakemake.config['args']['save_comments']):
        command.append("-C")
    
    command += ["-t", str(snakemake.threads), "-R", "@RG\\tID:"+snakemake.params.sample+"\\tSM:"+snakemake.params.sample, snakemake.input.ref, snakemake.input.fq1]

    if snakemake.config['in']['fq2'] != "None":
        command.append(snakemake.input.fq2)
    
    mccutils.run_command_stdout(command, snakemake.output[0], log=snakemake.log[0])
    print("<PROCESSING> read mapping complete")
        



if __name__ == "__main__":                
    main()