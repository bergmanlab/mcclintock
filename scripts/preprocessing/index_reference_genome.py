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
    try:
        log = snakemake.params.log
        print("<PROCESSING> making samtools and bwa index files for reference fasta")
        mccutils.run_command(["samtools", "faidx", snakemake.input.ref],log=log)
        mccutils.run_command(["bwa", "index", snakemake.input.ref], log=log)

        for out in snakemake.output:
            mccutils.check_file_exists(out)
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...unable to index (samtools, bwa) reference fasta, please check the formatting of:", snakemake.input.ref, file=sys.stderr)
        sys.exit(1) 


if __name__ == "__main__":                
    main()