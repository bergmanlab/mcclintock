import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    command = ["perl", snakemake.input.script, snakemake.input.ref_gff, snakemake.input.taxonomy, "Alias"]
    mccutils.run_command(command)
        

if __name__ == "__main__":                
    main()