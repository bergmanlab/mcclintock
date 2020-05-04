#!/usr/bin/env python3

import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.preprocessing.trimgalore as trimgalore


def main():
    fq1 = snakemake.input[0]
    fq2 = snakemake.params.fq2
    methods = snakemake.config['args']['methods'].split(",")
    processors = snakemake.config['args']['proc']
    mcc_out = snakemake.config['args']['out']
    run_id = snakemake.config['args']['run_id']
    log = snakemake.params.log

    print("<PROCESSING> prepping reads for McClintock....")
    # trims adaptors of input fastq(s)
    if "trimgalore" in methods:
        print("<PROCESSING> running trim_galore...log:"+log)
        if fq2 == "None":
            flags = trimgalore.SINGLE_END_FLAGS
            trimmedfq = run_trim_galore(fq1, run_id, log, mcc_out, cores=processors, flags=flags)
            mccutils.run_command(["touch", snakemake.output[1]])
        else:
            flags = trimgalore.PAIRED_END_FLAGS
            trimmedfq, trimmedfq2 = run_trim_galore(fq1, run_id, log, mcc_out, fq2=fq2, cores=processors, flags=flags)
            mccutils.run_command(["mv", trimmedfq2, snakemake.output[1]])
        
        mccutils.run_command(["mv", trimmedfq, snakemake.output[0]])
    
    # makes local unzipped copies of input fastq files
    else:
        if ".gz" in fq1:
            mccutils.run_command_stdout(["zcat",fq1], snakemake.output[0])
        else:
            mccutils.run_command(["cp", fq1, snakemake.output[0]])
        
        if fq2 == "None":
            mccutils.run_command(["touch", snakemake.output[1]])
        
        elif ".gz" in fq2:
            mccutils.run_command_stdout(["zcat",fq2], snakemake.output[1])
        
        else:
            mccutils.run_command(["cp", fq2, snakemake.output[1]])
    
    print("<PROCESSING> read setup complete")


def run_trim_galore(fq1, run_id, log, out, fq2=None, cores=1, flags=[]):
    mccutils.mkdir(out+"/results/")
    command = ['trim_galore'] + flags + ["-j", str(cores), "-o", out+"/results/trimgalore"]
    if fq2 is None:
        command.append(fq1)
    else:
        command += [fq1, fq2]
    
    mccutils.run_command(command, log=log)

    if fq2 is None:
        outfq = ""
        for f in os.listdir(out+"/results/trimgalore"):
            if "_trimmed.fq" in f:
                outfq = out+"/results/trimgalore/"+f
        if outfq != "":
            return outfq
        else:
            sys.exit("can't find trimgalore output fastq...exiting...\n")


    else:
        outfq1 = ""
        outfq2 = ""
        for f in os.listdir(out+"/results/trimgalore"):
            if "_val_1.fq" in f:
                outfq1 = out+"/results/trimgalore/"+f
            elif "_val_2.fq" in f:
                outfq2= out+"/results/trimgalore/"+f

        if outfq1 != "" and outfq2 != "":
            return outfq1, outfq2
        else:
            sys.exit("can't find trimgalore output fastq(s)...exiting...\n")

        


if __name__ == "__main__":                
    main()
