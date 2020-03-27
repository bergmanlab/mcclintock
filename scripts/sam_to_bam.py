import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    command = ["samtools","view", "-@", str(snakemake.threads), "-Sb", "-t", snakemake.input.ref_idx, snakemake.input.sam]
    mccutils.run_command_stdout(command, snakemake.output.tmp_bam)

    command = ["samtools", "sort", "-@", str(snakemake.threads), snakemake.output.tmp_bam, snakemake.output.bam.replace(".bam", "")]
    mccutils.run_command(command)

    command = ["samtools", "index", snakemake.output.bam]
    mccutils.run_command(command)

    command = ["samtools", "flagstat", snakemake.output.bam]
    mccutils.run_command_stdout(command, snakemake.output.flagstat)
        



if __name__ == "__main__":                
    main()