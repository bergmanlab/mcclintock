#!/usr/bin/env python3

import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    reference_fa = snakemake.input[0]
    consensus_TEs = snakemake.input[1]
    locations_gff = snakemake.input[2]
    family_tsv = snakemake.input[3]

    subprocess.call(["touch", snakemake.config['args']['out']+"/preprocessing.log"])

if __name__ == "__main__":                
    main()