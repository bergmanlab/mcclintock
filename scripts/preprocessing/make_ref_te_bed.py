import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    print("<PROCESSING> making reference TE bed file")
    with open(snakemake.input[0],"r") as i:
        with open(snakemake.output[0], "w") as o:
            for line in i:
                split_line = line.replace(";","\t").split("\t")
                name = split_line[9].split("=")[1]
                line = "\t".join([split_line[0], str(int(split_line[3])-1), split_line[4], name, ".", split_line[6]])
                o.write(line+"\n")
        

if __name__ == "__main__":                
    main()