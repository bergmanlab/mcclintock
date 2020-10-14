import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.jitterbug.jitterbug_run as config

def main():
    mccutils.log("jitterbug","Running jitterbug")

    reference_te_gff = snakemake.input.reference_tes
    bam = snakemake.input.bam

    out_dir = snakemake.params.out_dir
    script_dir = snakemake.params.script_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    threads = snakemake.threads

    out = snakemake.output.out

    out_gff, config = run_jitterbug(script_dir, bam, reference_te_gff, sample_name, out_dir, threads=threads, log=log)

    filter_jitterbug(script_dir, out_gff, config, sample_name, out, log=log)


def run_jitterbug(script_dir, bam, ref_te_gff, sample_name, out_dir, threads=1, log=None):
    command = [
        script_dir+"jitterbug.py", 
            "--mem", 
            "--numCPUs", str(threads), 
            "--output_prefix", out_dir+"/"+sample_name, 
            bam,
            ref_te_gff
    ]
    mccutils.run_command(command, log=log)

    return out_dir+"/"+sample_name+".TE_insertions_paired_clusters.gff3", out_dir+"/"+sample_name+".filter_config.txt"


def filter_jitterbug(script_dir, jitterbug_gff, filter_config, sample_name, filtered_gff, log=None):
    command = [
        script_dir+"tools/jitterbug_filter_results_func.py", 
            "-g", jitterbug_gff, 
            "-c", filter_config,
            "-o", filtered_gff
    ]

    mccutils.run_command(command, log=log)

if __name__ == "__main__":                
    main()