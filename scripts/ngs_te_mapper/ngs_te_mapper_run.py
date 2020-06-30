import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.ngs_te_mapper.ngs_te_mapper_run as config


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

    is_paired = True
    if snakemake.params.raw_fq2 == "None":
        is_paired = False
    
    command = ['Rscript', "--vanilla", script_dir+"/ngs_te_mapper.R", "genome="+reference_fasta, "teFile="+consensus_fasta, "tsd="+str(config.MAX_TSD), "thread="+str(threads), "output="+out_dir, "sourceCodeFolder="+script_dir]

    if is_paired:
        command.append("sample="+fastq1+";"+fastq2)
    else:
        command.append("sample="+fastq1)
    
    mccutils.log("ngs_te_mapper","running ngs_te_mapper", log=log)
    mccutils.run_command(command, log=log)
    mccutils.log("ngs_te_mapper","ngs_te_mapper run complete", log=log)


    raw_bed = ""
    for f in os.listdir(out_dir+"/bed_tsd/"):
        if "insertions.bed" in f:
            raw_bed = out_dir+"/bed_tsd/"+f

    mccutils.run_command(["cp", raw_bed, out_bed])

    mccutils.remove(out_dir+"/aligned_te/")

    mccutils.log("ngs_te_mapper","ngs_te_mapper run complete")

    



if __name__ == "__main__":                
    main()