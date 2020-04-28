import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils



def main():
    print("<POPOOLATIONTE> Running PopoolationTE...")
    ref_fasta = snakemake.input.ref_fasta
    taxonomy = snakemake.input.taxonomy
    te_gff = snakemake.input.te_gff
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log

    mccutils.run_command(["touch", snakemake.output[0]])

if __name__ == "__main__":                
    main()