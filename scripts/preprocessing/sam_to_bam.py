import os
import sys
import subprocess
import traceback
try:
    sys.path.append(snakemake.config['args']['mcc_path'])
    import scripts.mccutils as mccutils
except Exception as e:
    track = traceback.format_exc()
    print(track, file=sys.stderr)
    print("ERROR...unable to locate required external scripts at: "+snakemake.config['args']['mcc_path']+"/scripts/", file=sys.stderr)
    sys.exit(1)


def main():
    log = snakemake.params.log
    print("<PROCESSING> Converting sam to bam...log:"+log)

    try:
        command = ["samtools","view", "-@", str(snakemake.threads), "-Sb", "-t", snakemake.input.ref_idx, snakemake.input.sam]
        mccutils.run_command_stdout(command, snakemake.output.tmp_bam, log=log)
        mccutils.check_file_exists(snakemake.output.tmp_bam)

    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...unable convert sam to bam using SAMtools...sam file:", snakemake.input.sam, file=sys.stderr)
        sys.exit(1)


    try:
        command = ["samtools", "sort", "-@", str(snakemake.threads), snakemake.output.tmp_bam, snakemake.output.bam.replace(".bam", "")]
        mccutils.run_command(command, log=log)
        mccutils.check_file_exists(snakemake.output.bam)
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...falied to sort the bam file using samtools sort...bam file:", snakemake.output.tmp_bam, file=sys.stderr)
        sys.exit(1)

    try:
        command = ["samtools", "index", snakemake.output.bam]
        mccutils.run_command(command, log=log)

    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...falied to index the bam file using samtools index...bam file:", snakemake.output.bam, file=sys.stderr)
        sys.exit(1)


    try:
        command = ["samtools", "flagstat", snakemake.output.bam]
        mccutils.run_command_stdout(command, snakemake.output.flagstat, log=log)
        mccutils.check_file_exists(snakemake.output.flagstat)
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...falied to generate flagstat file using samtools flagstat...bam file:", snakemake.output.bam, file=sys.stderr)
        sys.exit(1)
        



if __name__ == "__main__":                
    main()