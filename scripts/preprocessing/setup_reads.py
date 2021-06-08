#!/usr/bin/env python3

import os
import sys
import subprocess
import traceback
from datetime import datetime
try:
    import importlib
    spec = importlib.util.spec_from_file_location("config", snakemake.params.config)
    trimgalore = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = trimgalore
    spec.loader.exec_module(trimgalore)
    sys.path.append(snakemake.config['args']['mcc_path'])
    import scripts.mccutils as mccutils
    from Bio import SeqIO
except Exception as e:
    track = traceback.format_exc()
    print(track, file=sys.stderr)
    print("ERROR...unable to locate required external scripts at: "+snakemake.config['args']['mcc_path']+"/scripts/ or "+snakemake.config['args']['mcc_path']+"/config/preprocessing/", file=sys.stderr)
    sys.exit(1)
    

class fileFormatError(Exception):
    def __init__(self, message):
        self.message = message
    pass

def main():
    fq1 = snakemake.input.fq1
    fq2 = snakemake.params.fq2
    methods = snakemake.params.methods.split(",")
    processors = snakemake.threads
    mcc_out = snakemake.params.out
    run_id = snakemake.params.run_id
    log = snakemake.params.log

    # now = datetime.now()
    # start = now.strftime("%Y-%m-%d %H:%M:%S")
    mccutils.log("processing", "prepping reads for McClintock")
    # trims adaptors of input fastq(s)
    trimmedfq = fq1
    trimmedfq2 = fq2

    try:
        check_fastqs(fq1, fq2, mcc_out, min_length=30, log=log)

        if "trimgalore" in methods:
            mccutils.log("processing", "running trim_galore", log=log)
            if fq2 == "None":
                trimmedfq = run_trim_galore(fq1, run_id, log, mcc_out, cores=processors, params=trimgalore.PARAMS["single_end"])
            else:
                trimmedfq, trimmedfq2 = run_trim_galore(fq1, run_id, log, mcc_out, fq2=fq2, cores=processors, params=trimgalore.PARAMS["paired_end"])
            
            run_multiqc(mcc_out+"/results/trimgalore/")
            
        
        # make unzipped copies in mcc input dir        
        make_copies(trimmedfq, trimmedfq2, snakemake.output[0], snakemake.output[1])
    
        # removes trimmed read files from trimgalore directory
        if trimmedfq != fq1:
            mccutils.remove(trimmedfq)
        if trimmedfq2 != fq2:
            mccutils.remove(trimmedfq2)

    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR processing of FastQ files failed...check that your FastQ files are formatted correctly...Exiting...", file=sys.stderr)
        mccutils.remove(snakemake.output[0])
        mccutils.remove(snakemake.output[1])
        sys.exit(1)


    # now = datetime.now()
    # end = now.strftime("%Y-%m-%d %H:%M:%S")
    # mccutils.log("setup_reads", "start: "+start)
    # mccutils.log("setup_reads", "end: "+end)

    mccutils.log("processing", "read setup complete")


def make_copies(fq1, fq2, fq1copy, fq2copy):
    if "gz" in fq1.split(".")[-1]:
        mccutils.run_command_stdout(["zcat",fq1], fq1copy)
    else:
        mccutils.run_command(["cp", fq1, fq1copy])
    
    if fq2 == "None":
        mccutils.run_command(["touch", fq2copy])
    
    elif "gz" in fq2.split(".")[-1]:
        mccutils.run_command_stdout(["zcat",fq2], fq2copy)
    
    else:
        mccutils.run_command(["cp", fq2, fq2copy])
    
    return fq1copy, fq2copy

def has_valid_read_lengths(fq1, fq2, min_length=30, paired=False):
    if paired:
        fqs_to_check = [fq1, fq2]
    else:
        fqs_to_check = [fq1]
    
    for x,fq in enumerate(fqs_to_check):
        has_valid_reads = False
        for record in SeqIO.parse(fq, "fastq"):
            if len(str(record.seq)) >= min_length:
                has_valid_reads = True
                break
        
        if not has_valid_reads:
            raise fileFormatError("fastq "+str(x+1)+" lacks any reads >= the minimum length of:"+str(min_length))

def has_valid_read_ids(fq1, fq2, log=None):
    passed = mccutils.run_command(["fastq_info", fq1, fq2], log=log, fatal=False)
    if not passed:
        raise fileFormatError("Paired fastq files failed validation, see: "+log+" for details")



def check_fastqs(fq1, fq2, out, min_length=30, log=None):
    mccutils.mkdir(out+"/tmp")
    if fq2 == "None":
        paired = False
    else:
        paired =True
    
    fq1, fq2 = make_copies(fq1, fq2, out+"/tmp/tmp_val_fq_1.fq", out+"/tmp/tmp_val_fq_2.fq")

    has_valid_read_lengths(fq1, fq2, min_length=min_length, paired=paired)

    if paired:
        has_valid_read_ids(fq1, fq2, log=log)


def run_trim_galore(fq1, run_id, log, out, fq2=None, cores=1, params={}):
    mccutils.mkdir(out+"/results/")
    command = ['trim_galore', "-j", str(cores), "-o", out+"/results/trimgalore"]

    for param in params.keys():
        if params[param] == True:
            command.append(param)

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

        file_exists = mccutils.check_file_exists(outfq)
        return outfq

    else:
        outfq1 = ""
        outfq2 = ""
        for f in os.listdir(out+"/results/trimgalore"):
            if "_val_1.fq" in f:
                outfq1 = out+"/results/trimgalore/"+f
            elif "_val_2.fq" in f:
                outfq2= out+"/results/trimgalore/"+f

        file_exists = mccutils.check_file_exists(outfq1)
        file_exists = mccutils.check_file_exists(outfq2)
        return outfq1, outfq2

def run_multiqc(trimgalore_dir):
    os.chdir(trimgalore_dir)
    mccutils.run_command(["multiqc","."])


if __name__ == "__main__":                
    main()
