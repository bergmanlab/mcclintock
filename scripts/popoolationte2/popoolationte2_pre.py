import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils



def main():
    mccutils.log("popoolationte2","setting up for PopoolationTE2")
    ref_fasta = snakemake.input.ref_fasta
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2
    log = snakemake.params.log
    out_dir = snakemake.params.out_dir
    threads = snakemake.threads

    # ensures intermediate files from previous runs are removed
    for f in os.listdir(out_dir):
        mccutils.remove(out_dir+"/"+f)

    index_fasta(ref_fasta, log=log)
    sam = map_reads(ref_fasta, fq1, fq2, out_dir, threads=threads, log=log)
    bam = sam_to_bam(sam, out_dir+"/tmp.bam", threads=threads, log=log)
    sorted_bam = sort_bam(bam, snakemake.output.bam, threads=threads, log=log)

    mccutils.remove(sam)
    mccutils.remove(bam)


def index_fasta(fasta, log=None):
    mccutils.log("popoolationte2","indexing reference fasta", log=log)
    command = ["bwa", "index", fasta]
    mccutils.run_command(command, log=log)

def map_reads(ref, fq1, fq2, out, threads=1, log=None):
    mccutils.log("popoolationte2","mapping reads", log=log)
    sam = out+"/"+"mapped.sam"
    mccutils.run_command_stdout(["bwa","mem", "-M", "-t", str(threads), ref, fq1, fq2], sam, log=log)
    return sam

def sam_to_bam(sam, bam, threads=1, log=None):
    mccutils.log("popoolationte2","converting SAM to BAM", log=log)
    mccutils.run_command_stdout(["samtools", "view","-@",str(threads), "-Sb", sam], bam, log=log)
    return bam

def sort_bam(bam, sorted_bam, threads=1, log=None):
    mccutils.log("popoolationte2","sorting BAM", log=log)
    mccutils.run_command(["samtools","sort","-@", str(threads), bam, "-o", sorted_bam], log=log)
    return sorted_bam


if __name__ == "__main__":                
    main()