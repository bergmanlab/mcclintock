import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    mccutils.log("processing","creating ngs_te_mapper2 reference TE gff")
    in_gff = snakemake.input.locations
    taxonomy = snakemake.input.taxonomy
    out_gff = snakemake.output[0]

    te_families = {}
    with open(taxonomy,"r") as tax:
        for line in tax:
            split_line = line.split("\t")
            te_families[split_line[0]] = split_line[1].replace("\n","")
    
    with open(in_gff, "r") as ingff:
        with open(out_gff,"w") as outgff:
            for line in ingff:
                split_line = line.split("\t")
                split_line[8] = "Target="+te_families[split_line[2]]
                line = "\t".join(split_line)
                outgff.write(line+"\n")
    
    mccutils.log("processing","ngs_te_mapper2 reference TE gff created")

if __name__ == "__main__":                
    main()