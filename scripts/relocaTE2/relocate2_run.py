import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.relocate2.relocate2_run as config


def main():
    reference = snakemake.input.reference
    te_seqs = snakemake.input.te_seqs
    rm_out = snakemake.input.rm_out
    fq1 = snakemake.input.fastq1
    fq2 = snakemake.input.fastq2
    bam = snakemake.input.bam
    log = snakemake.params.log
    raw_fq2 = snakemake.params.raw_fq2
    is_paired = True
    if snakemake.params.raw_fq2 == "None":
        is_paired = False
    
    with open(log,"a") as l:
        l.write("reference fasta: "+reference+"\n")
        l.write("TE fasta: "+te_seqs+"\n")
        l.write("repeatmasker out: "+rm_out+"\n")
        l.write("BAM: "+bam+"\n")
        l.write("fastq1: "+fq1+"\n")
        if is_paired:
            l.write("fastq2: "+fq2+"\n")


    threads = snakemake.threads
    out_dir = snakemake.params.out_dir
    median_insert_size_file = snakemake.input.median_insert_size

    median_insert_size = get_median_insert_size(median_insert_size_file)
    fq_dir = os.path.dirname(fq1)
    output = subprocess.Popen(["which", "relocaTE2.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    script = output.stdout.read()
    script = script.decode()
    script = script.replace("\n","")




    print("<RELOCATE> Running RelocaTE2...")
    command =  [
        "python2", script, 
                        "-t", te_seqs,
                        "-g", reference,
                        "-r", rm_out,
                        "-o", out_dir,
                        "-s", str(median_insert_size),
                        "--run",
                        "-v", "4",
                        "-c", str(threads),
                        "--aligner", config.RELOCATE2["aligner"],
                        "--len_cut_match", str(config.RELOCATE2["len_cut_match"]),
                        "--len_cut_trim", str(config.RELOCATE2["len_cut_trim"]),
                        "--mismatch", str(config.RELOCATE2["mismatch"]),
                        "--mismatch_junction", str(config.RELOCATE2["mismatch_junction"])
    ]

    if is_paired:
        command += ["-1", "_1", "-2", "_2", "-d", fq_dir+"/"]
    
    else:
        mccutils.mkdir(out_dir+"/tmp")
        mccutils.mkdir(out_dir+"/tmp/fastq")
        mccutils.run_command(["cp", fq1, out_dir+"/tmp/fastq"])
        command += ["-u", "_1", "-d", out_dir+"/tmp/fastq"]


    print(" ".join(command))
    mccutils.run_command(command, log=log)


def get_median_insert_size(infile):
    median_insert_size = 0
    with open(infile,"r") as inf:
        for line in inf:
            insert = line.split("=")[1]
            insert = insert.replace("\n","")
            median_insert_size = int(float(insert))
    
    return median_insert_size

if __name__ == "__main__":                
    main()