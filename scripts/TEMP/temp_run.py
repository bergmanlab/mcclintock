import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.TEMP.temp_run as config



def main():
    print("<TEMP> Running TEMP Module...")
    bam = snakemake.input.bam
    twobit = snakemake.input.twobit
    consensus = snakemake.input.consensus
    ref_te_bed = snakemake.input.ref_te_bed
    taxonomy = snakemake.input.taxonomy
    median_insert_size_file = snakemake.input.median_insert_size
    threads = snakemake.threads
    out_dir = snakemake.params.out_dir
    scripts_dir = snakemake.params.scripts_dir
    log = snakemake.params.log

    sample_name = snakemake.params.sample_name

    median_insert_size = get_median_insert_size(median_insert_size_file)

    run_temp_insertion(bam, scripts_dir, consensus, ref_te_bed, taxonomy, median_insert_size, threads, out_dir, log)

    run_temp_absence(bam, scripts_dir, consensus, ref_te_bed, twobit, taxonomy, median_insert_size, threads, out_dir, log)



def get_median_insert_size(infile):
    median_insert_size = 0
    with open(infile,"r") as inf:
        for line in inf:
            insert = line.split("=")[1]
            insert = insert.replace("\n","")
            median_insert_size = int(float(insert))
    
    return median_insert_size

    
def run_temp_insertion(bam, scripts, consensus, te_bed, taxonomy, median_insert_size, threads, out, log):
    print("<TEMP> Running TEMP non-reference insertion prediction...")
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


def run_temp_absence(bam, scripts, consensus, te_bed, twobit, taxonomy, median_insert_size, threads, out, log):
    print("<TEMP> Running TEMP reference insertion absence prediction...")
    command = [
        "bash", scripts+"TEMP_Absence.sh",
            "-i", bam,
            "-s", scripts,
            "-r", te_bed,
            "-t", twobit,
            "-f", str(median_insert_size),
            "-c", str(threads),
            "-o", out
    ]

    mccutils.run_command(command, log=log)


if __name__ == "__main__":                
    main()