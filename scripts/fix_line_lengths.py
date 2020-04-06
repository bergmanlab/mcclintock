#!/usr/bin/env python3

import sys
from Bio import SeqIO
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.fix_fasta as fix_fasta
import modules.mccutils as mccutils


def main():
    print("<PROCESSING> running rule: fix_line_lengths")
    length = 80

    fasta1 = snakemake.input[0]
    lines = fix_fasta.fix_fasta_lines(fasta1, length)
    write_fasta(lines, snakemake.output[0])

    fasta2 = snakemake.input[1]
    lines = fix_fasta.fix_fasta_lines(fasta2, length)
    write_fasta(lines, snakemake.output[1])

    if snakemake.params.coverage_fasta == "None":
        mccutils.run_command(["touch",snakemake.output[2]])
    else:
        fasta3 = snakemake.params.coverage_fasta
        lines = fix_fasta.fix_fasta_lines(fasta3, length)
        write_fasta(lines, snakemake.output[2])
    
    print("<PROCESSING> rule: fix_line_lengths Complete")
        


def write_fasta(lines, fasta):
    with open(fasta, "w") as out:
        for line in lines:
            out.write(line+"\n")

if __name__ == "__main__":                
    main()