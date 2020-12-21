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
    jar = snakemake.params.jar
    log = snakemake.params.log
    out_dir = snakemake.params.out_dir
    threads = snakemake.threads

    # ensures intermediate files from previous runs are removed
    for f in os.listdir(out_dir):
        mccutils.remove(out_dir+"/"+f)

    mccutils.mkdir(out_dir+"/tmp")
    index_fasta(ref_fasta, log=log)
    fq1 = format_fastq(fq1, out_dir+"/reads_1.fastq", log=log)
    fq2 = format_fastq(fq2, out_dir+"/reads_2.fastq", log=log)
    sam1 = map_reads(ref_fasta, fq1, out_dir+"/mapped_1.sam", threads=threads, log=log)
    sam2 = map_reads(ref_fasta, fq2, out_dir+"/mapped_2.sam", threads=threads, log=log)
    bam = sam_to_bam(jar, fq1, fq2, sam1, sam2, snakemake.output.bam, out_dir, threads=threads, log=log)
    mccutils.remove(out_dir+"/tmp")


def index_fasta(fasta, log=None):
    mccutils.log("popoolationte2","indexing reference fasta", log=log)
    command = ["bwa", "index", fasta]
    mccutils.run_command(command, log=log)

def format_fastq(fq, out_fq, log=None):
    mccutils.log("popoolationte2","formatting fastq read names", log=log)
    with open(fq,"r") as inf:
        with open(out_fq,"w") as of:
            for ln,line in enumerate(inf, start=4):
                if (ln%4) == 0:
                    line = line.replace("\n","")
                    line = line.split(" ")[0]
                    line += "\n"
                of.write(line)
    
    return out_fq


def map_reads(ref, fq, outsam, threads=1, log=None):
    mccutils.log("popoolationte2","mapping reads", log=log)
    mccutils.run_command_stdout(["bwa","bwasw", "-t", str(threads), ref, fq], outsam, log=log)
    return outsam

def sam_to_bam(jar, fq1, fq2, sam1, sam2, bam, out_dir, threads=1, log=None):
    mccutils.log("popoolationte2","converting SAM to BAM", log=log)
    mccutils.run_command(["java", "-Djava.io.tmpdir="+out_dir+"/tmp", "-jar",jar, "se2pe", 
                                              "--fastq1", fq1,
                                              "--fastq2", fq2,
                                              "--bam1", sam1,
                                              "--bam2", sam2,
                                              "--sort",
                                              "--output", bam], log=log)
    return bam



if __name__ == "__main__":                
    main()