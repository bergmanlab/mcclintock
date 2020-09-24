#!/usr/bin/env python3

import os
import sys
import subprocess
import traceback
try:
    sys.path.append(snakemake.config['args']['mcc_path'])
    import scripts.mccutils as mccutils
    import config.preprocessing.trimgalore as trimgalore
except Exception as e:
    track = traceback.format_exc()
    print(track, file=sys.stderr)
    print("ERROR...unable to locate required external scripts at: "+snakemake.config['args']['mcc_path']+"/scripts/ or "+snakemake.config['args']['mcc_path']+"/config/preprocessing/", file=sys.stderr)
    sys.exit(1)
    


def main():
    fq1 = snakemake.input.fq1
    fq2 = snakemake.params.fq2
    methods = snakemake.params.methods.split(",")
    processors = snakemake.threads
    mcc_out = snakemake.params.out
    run_id = snakemake.params.run_id
    log = snakemake.params.log

    mccutils.log("processing", "prepping reads for McClintock")
    try:
        # trims adaptors of input fastq(s)
        trimmedfq = fq1
        trimmedfq2 = fq2
        if "trimgalore" in methods:
            mccutils.log("processing", "running trim_galore", log=log)
            if fq2 == "None":
                flags = trimgalore.SINGLE_END_FLAGS
                trimmedfq = run_trim_galore(fq1, run_id, log, mcc_out, cores=processors, flags=flags)
            else:
                flags = trimgalore.PAIRED_END_FLAGS
                trimmedfq, trimmedfq2 = run_trim_galore(fq1, run_id, log, mcc_out, fq2=fq2, cores=processors, flags=flags)
            
            run_multiqc(mcc_out+"/results/trimgalore/")
            
        
        # make unzipped copies in mcc input dir        
        if "gz" in trimmedfq.split(".")[-1]:
            mccutils.run_command_stdout(["zcat",trimmedfq], snakemake.output[0])
        else:
            mccutils.run_command(["cp", trimmedfq, snakemake.output[0]])
        
        if trimmedfq2 == "None":
            mccutils.run_command(["touch", snakemake.output[1]])
        
        elif "gz" in trimmedfq2.split(".")[-1]:
            mccutils.run_command_stdout(["zcat",trimmedfq2], snakemake.output[1])
        
        else:
            mccutils.run_command(["cp", trimmedfq2, snakemake.output[1]])
    
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

    mccutils.log("processing", "read setup complete")


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
