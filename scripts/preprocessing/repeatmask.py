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

    run_repeatmasker(reference, ref_name, te_seqs, threads, log, output, out_dir)





def run_repeatmasker(reference, ref_name, te_seqs, threads, log, outfile, outdir):
    tmp_dir = outdir+"/tmp/repeatmasker"
    mccutils.mkdir(tmp_dir)
    os.chdir(tmp_dir)

    command = ["RepeatMasker","-pa", str(threads), "-lib", te_seqs, "-dir", tmp_dir, "-s", "-gff", "-nolow", "-no_is", reference]
    mccutils.run_command(command, log=log)

    os.chdir(outdir)

    rm_out = ""
    for f in os.listdir(tmp_dir):
        if ".out" in f:
            rm_out = tmp_dir+"/"+f

    mccutils.run_command(["mv", rm_out, outfile])




if __name__ == "__main__":                
    main()