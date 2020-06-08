import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils



def main():
    print("<POPOOLATIONTE2> Setting up for PopoolationTE2...")
    ref_fasta = snakemake.input.ref_fasta
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2
    log = snakemake.params.log
    out_dir = snakemake.params.out_dir
    threads = snakemake.threads

    index_fasta(ref_fasta, log=log)
    sam = map_reads(ref_fasta, fq1, fq2, out_dir, threads=threads, log=log)
    bam = sam_to_bam(sam, out_dir+"/tmp.bam", threads=threads, log=log)
    sorted_bam = sort_bam(bam, snakemake.output.bam, threads=threads, log=log)


def index_fasta(fasta, log=None):
    print("<POPOOLATIONTE2> Indexing reference fasta...")
    command = ["bwa", "index", fasta]
    mccutils.run_command(command, log=log)

def map_reads(ref, fq1, fq2, out, threads=1, log=None):
    print("<POPOOLATIONTE2> Mapping reads...")
    sam = out+"/"+"mapped.sam"
    mccutils.run_command_stdout(["bwa","bwasw", "-t", str(threads), ref, fq1, fq2], sam, log=log)
    return sam

def sam_to_bam(sam, bam, threads=1, log=None):
    print("<POPOOLATIONTE2> Converting SAM to BAM...")
    mccutils.run_command_stdout(["samtools", "view","-@",str(threads), "-Sb", sam], bam, log=log)
    return bam

def sort_bam(bam, sorted_bam, threads=1, log=None):
    print("<POPOOLATIONTE2> Sorting BAM...")
    mccutils.run_command(["samtools","sort","-@", str(threads), "-f", bam, sorted_bam], log=log)
    return sorted_bam


if __name__ == "__main__":                
    main()