import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    print("<PROCESSING> formatting the name of consensus TE fasta headers for compatibility with relocaTE...")
    with open(snakemake.input[0],"r") as infa:
        with open(snakemake.output[0], "w") as outfa:
            for line in infa:
                if ">" in line:
                    line = line.replace("\n","")
                    line += " TSD=UNK\n"

                outfa.write(line)
        

if __name__ == "__main__":                
    main()