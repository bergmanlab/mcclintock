import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    mccutils.log("processing","creating relocaTE reference TE gff")
    taxonomy = {}
    with open(snakemake.input[1],"r") as tax:
        for line in tax:
            split_line = line.split("\t")
            taxonomy[split_line[0]] = split_line[1].replace("\n","")
    
    with open(snakemake.input[0], "r") as ingff:
        with open(snakemake.output[0],"w") as outgff:
            for line in ingff:
                split_line = line.split("\t")
                split_line[2] = taxonomy[split_line[2]]
                line = "\t".join(split_line)
                outgff.write(line)
    
    mccutils.log("processing","relocaTE reference TE gff created")

if __name__ == "__main__":                
    main()