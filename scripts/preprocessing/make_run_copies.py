import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    print("<PROCESSING> making run copies of input files")
    te_locations = snakemake.config['in']['locations']
    taxonomy = snakemake.config['in']['taxonomy']

    if te_locations == "None":
        mccutils.run_command(["touch", snakemake.output[0]])
    else:
        mccutils.run_command(["cp",te_locations, snakemake.output[0]])
    
    if taxonomy == "None":
        mccutils.run_command(["touch", snakemake.output[1]])
    else:
        mccutils.run_command(["cp", taxonomy, snakemake.output[1]])


if __name__ == "__main__":                
    main()