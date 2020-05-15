#!/usr/bin/env python3

import sys
from Bio import SeqIO
import traceback
try:
    sys.path.append(snakemake.config['args']['mcc_path'])
    import scripts.fix_fasta as fix_fasta
    import scripts.mccutils as mccutils
except Exception as e:
    track = traceback.format_exc()
    print(track, file=sys.stderr)
    print("ERROR...unable to locate required external scripts at: "+snakemake.config['args']['mcc_path']+"/scripts/", file=sys.stderr)
    sys.exit(1)


def main():
    print("<PROCESSING> running rule: fix_line_lengths")
    fastas = []
    try:
        length = 80
        fasta1 = snakemake.input[0]
        lines = fix_fasta.fix_fasta_lines(fasta1, length)
        write_fasta(lines, snakemake.output[0])

        fasta2 = snakemake.input[1]
        lines = fix_fasta.fix_fasta_lines(fasta2, length)
        write_fasta(lines, snakemake.output[1])

        fastas = [fasta1, fasta2]
        if snakemake.params.coverage_fasta == "None":
            mccutils.run_command(["touch",snakemake.output[2]])
        else:
            fasta3 = snakemake.params.coverage_fasta
            fastas.append(fasta3)
            lines = fix_fasta.fix_fasta_lines(fasta3, length)
            write_fasta(lines, snakemake.output[2])
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...failed to fix the line lengths of the input fasta files, check the formatting of :"+",".join(fastas), file=sys.stderr)
        mccutils.remove(snakemake.output[0])
        mccutils.remove(snakemake.output[1])
        mccutils.remove(snakemake.output[2])
        sys.exit(1)        

    print("<PROCESSING> rule: fix_line_lengths Complete")
        


def write_fasta(lines, fasta):
    with open(fasta, "w") as out:
        for line in lines:
            out.write(line+"\n")

if __name__ == "__main__":                
    main()