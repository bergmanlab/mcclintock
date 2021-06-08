#!/usr/bin/env python3

import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    reference = snakemake.input.reference
    te_seqs = snakemake.input.te_seqs

    ref_name = snakemake.params.ref_name
    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    threads = snakemake.threads

    output = snakemake.output.rm_out

    mccutils.log("processing","Running RepeatMasker", log=log)
    run_repeatmasker(reference, ref_name, te_seqs, threads, log, output, out_dir)
    mccutils.log("processing","Repeatmasker complete")





def run_repeatmasker(reference, ref_name, te_seqs, threads, log, outfile, outdir):
    tmp_dir = outdir+"/tmp/repeatmasker"
    mccutils.remove(tmp_dir)
    mccutils.mkdir(tmp_dir)
    os.chdir(tmp_dir)

    command = ["RepeatMasker","-pa", str(threads), "-lib", te_seqs, "-dir", tmp_dir, "-s", "-nolow", "-no_is", reference]
    mccutils.run_command(command, log=log)

    os.chdir(outdir)

    rm_out = ""
    for f in os.listdir(tmp_dir):
        if "fasta.out" in f and f[-9:] == "fasta.out":
            rm_out = tmp_dir+"/"+f

    if rm_out == "":
        sys.exit("can't find Repeatmasker output in:"+tmp_dir+"\n")

    mccutils.run_command(["mv", rm_out, outfile])




if __name__ == "__main__":                
    main()