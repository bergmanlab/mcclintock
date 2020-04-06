import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils
import config.temp as config



def main():
    print("<TEMP> Running TEMP Module...")
    bam = snakemake.input.bam
    twobit = snakemake.input.twobit
    consensus = snakemake.input.consensus
    ref_te_bed = snakemake.input.ref_te_bed
    taxonomy = snakemake.input.taxonomy
    median_insert_size_file = snakemake.input.median_insert_size
    te_gff = snakemake.input.te_gff
    sample_name = snakemake.params.sample
    threads = snakemake.threads
    out_dir = snakemake.params.out_dir
    scripts_dir = snakemake.params.scripts_dir
    log = snakemake.params.log

    median_insert_size = get_median_insert_size(median_insert_size_file)

    run_temp_insertion(bam, scripts_dir, consensus, ref_te_bed, taxonomy, median_insert_size, threads, out_dir, log)

    mccutils.run_command(["touch", out_dir+"TEMP.log"])


def get_median_insert_size(infile):
    median_insert_size = 0
    with open(infile,"r") as inf:
        for line in inf:
            insert = line.split("=")[1]
            insert = insert.replace("\n","")
            median_insert_size = int(insert)
    
    return median_insert_size

    
def run_temp_insertion(bam, scripts, consensus, te_bed, taxonomy, median_insert_size, threads, out, log):
    command = [
        "bash", scripts+"TEMP_Insertion.sh", 
            "-x", str(config.TEMP_Insertion['x']), 
            "-i", bam, 
            "-s", scripts, 
            "-r", consensus, 
            "-t", te_bed, 
            "-u", taxonomy, 
            "-m", str(config.TEMP_Insertion['m']),
            "-f", str(median_insert_size),
            "-c", str(threads),
            "-o", out
    ]

    mccutils.run_command(command, log=log)



if __name__ == "__main__":                
    main()