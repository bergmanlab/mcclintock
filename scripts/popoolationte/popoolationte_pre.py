import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils



def main():
    print("<POPOOLATIONTE PREPROCESSING> Running PopoolationTE preprocessing steps...")
    ref_fasta = snakemake.input.ref_fasta
    taxonomy = snakemake.input.taxonomy
    te_gff = snakemake.input.te_gff
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log

    threads = snakemake.threads

    fq1,fq2 = format_read_names(fq1, fq2, sample_name, out_dir)

    index_fasta(ref_fasta, log=log)
    map_reads(fq1, ref_fasta, snakemake.output[0], threads=threads, log=log)
    map_reads(fq2, ref_fasta, snakemake.output[1], threads=threads, log=log)

    print("<POPOOLATIONTE PREPROCESSING> PopoolationTE preprocessing complete")

def format_read_names(fq1, fq2, sample_name, out_dir):
    outfq1 = out_dir+"/reads1.fastq"
    outfq2 = out_dir+"/reads2.fastq"

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

def map_reads(fq, fasta, outfile, threads=1, log=None):
    command = ["bwa", "bwasw", "-t", str(threads), fasta, fq]
    mccutils.run_command_stdout(command, outfile, log=log)

def get_read_length(fq1, fq2):
    read1_length = 0
    with open(fq1, "r") as fq:
        for l, line in fq:
            if l == 1:
                read1_length = len(line.replace("\n",""))
            elif l > 1:
                break
    
    read2_length = 0
    with open(fq2, "r") as fq:
        for l, line in fq:
            if l == 1:
                read2_length = len(line.replace("\n",""))
            elif l > 1:
                break
    
    read_length = int((read1_length + read2_length)//2)

    return(read_length)


if __name__ == "__main__":                
    main()