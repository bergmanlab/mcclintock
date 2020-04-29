import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    print("<PROCESSING> adding fake chromosomes if chrom # < 5, required by TE-locate")
    chromosomes = []
    with open(snakemake.input[0],"r") as infa:
        with open(snakemake.output[0],"w") as outfa:
            for line in infa:
                outfa.write(line)
                if ">" in line:
                    chromosomes.append(line)

    if len(chromosomes) < 5:
        diff = 5 - len(chromosomes)        
        with open(snakemake.output[0],"a") as out: 
            for i in range(1,diff+1):
                out.write(">fixforTElocate"+str(i)+"\n")
                out.write("ACGT\n")
        

if __name__ == "__main__":                
    main()