import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils



def main():
    print("<POPOOLATIONTE PREPROCESSING> Running PopoolationTE preprocessing steps...")
    ref_fasta = snakemake.input.ref_fasta
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    script_dir = snakemake.params.script_dir

    threads = snakemake.threads

    print("\tformatting read names...")
    fq1,fq2 = format_read_names(fq1, fq2, sample_name, out_dir)
    print("\tindexing popoolationTE reference fasta...")
    index_fasta(ref_fasta, log=log)
    print("\tmapping reads from "+fq1)
    sam1 = map_reads(fq1, ref_fasta, threads=threads, log=log)
    print("\tmapping reads from "+fq2)
    sam2 = map_reads(fq2, ref_fasta, threads=threads, log=log)
    print("\tcombining alignments...")
    combined_sam = combine_alignments(sam1, sam2, fq1, fq2, script_dir, out_dir, log=log)
    print("\tsorting sam file...")
    bam = sam_to_bam(combined_sam, threads=threads, log=log)
    sorted_bam = sort_bam(bam, threads=threads, log=log)
    sorted_sam = bam_to_sam(sorted_bam, threads=threads, log=log)

    files_to_remove = [sam1, sam2, combined_sam, bam, sorted_bam]

    for f in files_to_remove:
        mccutils.remove(f)

    print("<POPOOLATIONTE PREPROCESSING> PopoolationTE preprocessing complete")

def format_read_names(fq1, fq2, sample_name, out_dir):
    outfq1 = out_dir+"reads1.fastq"
    outfq2 = out_dir+"reads2.fastq"

    with open(outfq1,"w") as outfq:
        with open(fq1,"r") as infq:
            read_num = 1
            for l,line in enumerate(infq):
                if l%4 == 0:
                    outfq.write("@"+sample_name+str(read_num)+"/1\n")
                    read_num += 1
                else:
                    outfq.write(line)
    
    with open(outfq2,"w") as outfq:
        with open(fq2,"r") as infq:
            read_num = 1
            for l,line in enumerate(infq):
                if l%4 == 0:
                    outfq.write("@"+sample_name+str(read_num)+"/2\n")
                    read_num += 1
                else:
                    outfq.write(line)

    return outfq1, outfq2

def index_fasta(fasta, log=None):
    command = ["bwa", "index", fasta]
    mccutils.run_command(command, log=log)

def map_reads(fq, fasta, threads=1, log=None):
    outfile = fq.split(".")
    outfile[-1] = "sam"
    outfile = ".".join(outfile)

    command = ["bwa", "bwasw", "-t", str(threads), fasta, fq]
    mccutils.run_command_stdout(command, outfile, log=log)

    return outfile


def combine_alignments(sam1, sam2, fq1, fq2, script_path, out, log=None):
    out_sam = out+"combined.sam"
    command = ["perl", script_path+"samro.pl", "--sam1", sam1, "--sam2", sam2, "--fq1", fq1, "--fq2", fq2, "--output", out_sam]
    mccutils.run_command(command, log=log)
    return out_sam


def sam_to_bam(sam, threads=1, log=None):
    bam = sam.split(".")
    bam[-1] = "bam"
    bam = ".".join(bam)

    command = ["samtools","view", "-Sb", "-@", str(threads), sam]
    mccutils.run_command_stdout(command, bam, log=log)

    return bam

def sort_bam(bam, threads=1, log=None):
    sorted_bam = bam.split(".")
    sorted_bam[-1] = "sorted.bam"
    sorted_bam = ".".join(sorted_bam)
    command = ["samtools", "sort", bam, "-@", str(threads), "-o", sorted_bam]
    mccutils.run_command_stdout(command, sorted_bam, log=log)

    return sorted_bam

def bam_to_sam(bam, threads=1, log=None):
    sam = bam.split(".")
    sam[-1] = "sam"
    sam = ".".join(sam)

    command = ["samtools", "view", "-@", str(threads), bam]
    mccutils.run_command_stdout(command, sam, log=log)

    return sam


if __name__ == "__main__":                
    main()