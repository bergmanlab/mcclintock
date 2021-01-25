import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.ngs_te_mapper2.ngs_te_mapper2_run as config


def main():
    consensus_fasta = snakemake.input.consensus_fasta
    reference_fasta = snakemake.input.reference_fasta
    fastq1 = snakemake.input.fastq1
    fastq2 = snakemake.input.fastq2

    log = snakemake.params.log
    with open(log,"a") as l:
        l.write("consensus fasta: "+consensus_fasta+"\n")
        l.write("reference fasta: "+reference_fasta+"\n")
        l.write("fastq1: "+fastq1+"\n")
        l.write("fastq2: "+fastq2+"\n")


    threads = snakemake.threads
    sample_name = snakemake.params.sample_name
    script_dir = snakemake.params.script_dir
    out_dir = snakemake.params.out_dir
    out_bed = snakemake.output[0]

    # ensures intermediate files from previous runs are removed
    for f in os.listdir(out_dir):
        mccutils.remove(out_dir+"/"+f)

    is_paired = True
    if snakemake.params.raw_fq2 == "None":
        is_paired = False
    
    command = [
        'python', script_dir+"/ngs_te_mapper.py", 
            "-r", reference_fasta, 
            "-l", consensus_fasta, 
            "--tsd_max", str(config.MAX_TSD), 
            "-t", str(threads), 
            "-o", out_dir, 
            "--keep_files",
            "-p", sample_name,
            "--af",
            "--min_af", str(config.MIN_ALLELE_FREQUENCY),
            "-f"
    ]

    if is_paired:
        command.append(fastq1+","+fastq2)
    else:
        command.append(fastq1)
    
    mccutils.log("ngs_te_mapper2","running ngs_te_mapper2", log=log)
    mccutils.run_command(command, log=log)
    mccutils.log("ngs_te_mapper2","ngs_te_mapper2 run complete", log=log)


    mccutils.log("ngs_te_mapper2","ngs_te_mapper2 run complete")

    



if __name__ == "__main__":                
    main()