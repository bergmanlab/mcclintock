import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.teflon.teflon_run as config

def main():
    mccutils.log("teflon","Running TEFLoN")

    consensus = snakemake.input.consensus
    reference_genome = snakemake.input.reference_genome
    ref_bed = snakemake.input.ref_bed
    teflon_taxonomy = snakemake.input.teflon_taxonomy
    bam = snakemake.input.bam

    threads = snakemake.threads
    out_dir = snakemake.params.out_dir
    script_dir = snakemake.params.script_dir
    log = snakemake.params.log

    sample_table = make_sample_table(out_dir, bam)

    run_teflon(
        script_dir, 
        out_dir, 
        sample_table, 
        threads=threads, 
        log=log, 
        quality_threshold=config.PARAMETERS['q'],
        stdev=config.PARAMETERS['sd'],
        cov=config.PARAMETERS['cov'],
        te_support1=config.PARAMETERS['n1'],
        te_support2=config.PARAMETERS['n2'],
        read_count_lower_threshold=config.PARAMETERS['lt'],
        read_count_higher_threshold=config.PARAMETERS['ht']
    )


def make_sample_table(out_dir, bam):
    with open(out_dir+"samples.tsv", "w") as samples:
        samples.write(bam+"\tsample\n")

    return out_dir+"samples.tsv"

def run_teflon(script_dir, out_dir, sample_file, threads=1, log=None, quality_threshold=20, stdev=None, cov=None, te_support1=1, te_support2=1, read_count_lower_threshold=1, read_count_higher_threshold=None):

    command = [
        "python", script_dir+"teflon.v0.4.py",
        "-wd", out_dir,
        "-d", out_dir+"teflon.prep_TF/",
        "-s", sample_file,
        "-i", "sample",
        "-l1", 'family',
        "-l2", 'family',
        "-t", str(threads),
        "-q", str(quality_threshold)
    ]

    if stdev is not None:
        command += ["-sd", str(stdev)]
    
    if cov is not None:
        command += ["-cov", str(cov)]

    mccutils.run_command(command, log=log)


    command = [
        "python", script_dir+"teflon_collapse.py",
        "-wd", out_dir,
        "-d", out_dir+"teflon.prep_TF/",
        "-s", sample_file,
        "-t", str(threads),
        "-n1", str(te_support1),
        "-n2", str(te_support2),
        "-q", str(quality_threshold)
    ]

    mccutils.run_command(command, log=log)

    command = [
        "python", script_dir+"teflon_count.py",
        "-wd", out_dir,
        "-d", out_dir+"teflon.prep_TF/",
        "-s", sample_file,
        "-i", "sample",
        "-l2", "family",
        "-t", str(threads),
        "-q", str(quality_threshold)
    ]

    mccutils.run_command(command, log=log)

    command = [
        "python", script_dir+"teflon_genotype.py",
        "-wd", out_dir,
        "-d", out_dir+"teflon.prep_TF/",
        "-s", sample_file,
        "-lt", str(read_count_lower_threshold),
        "-dt", "pooled"
    ]

    if read_count_higher_threshold is not None:
        command += ["-ht", str(read_count_higher_threshold)]

    mccutils.run_command(command, log=log)

if __name__ == "__main__":                
    main()