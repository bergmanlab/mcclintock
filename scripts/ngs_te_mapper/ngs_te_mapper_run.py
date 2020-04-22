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
    threads = snakemake.threads
    log = snakemake.params.log
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
    
    print("<NGS_TE_MAPPER> Running ngs_te_mapper...")
    mccutils.run_command(command, log=log)

    fq1_basename = mccutils.get_base_name(fastq1)
    fq2_basename = mccutils.get_base_name(fastq2)
    if is_paired:
        raw_bed = out_dir+"/bed_tsd/"+fq1_basename+"_"+fq2_basename+"insertions.bed"
    else:
        raw_bed = out_dir+"/bed_tsd/"+fq1_basename+"insertions.bed"

    mccutils.run_command(["cp", raw_bed, out_bed])

    print("<NGS_TE_MAPPER> ngs_te_mapper run complete...")

    



if __name__ == "__main__":                
    main()