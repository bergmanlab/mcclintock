import os
import sys
import subprocess
from Bio import SeqIO
import traceback
try:
    sys.path.append(snakemake.config['args']['mcc_path'])
    import scripts.mccutils as mccutils
    import scripts.fix_fasta as fix_fasta
except Exception as e:
    track = traceback.format_exc()
    print(track, file=sys.stderr)
    print("ERROR...unable to locate required external scripts at: "+snakemake.config['args']['mcc_path']+"/scripts/", file=sys.stderr)
    sys.exit(1)


def main():
    ref_fasta = snakemake.input.ref_fasta
    te_gff = snakemake.input.te_gff
    mcc_out = snakemake.params.mcc_out
    run_id = snakemake.params.run_id
    log = snakemake.params.log
    augment = snakemake.params.augment
    chromosomes = snakemake.params.chromosomes.split(",")
    ref_te_fasta = snakemake.output.ref_te_fasta


    ref_tes = get_ref_te_fasta(ref_fasta, te_gff, run_id, log, mcc_out)

    mccutils.run_command(["mv", ref_tes, ref_te_fasta])


def get_ref_te_fasta(reference, ref_tes_gff, run_id, log, out):
    try:
        ref_te_fasta = out+"/tmp/"+run_id+"tmpreferencetes.fasta"
        command = ["bedtools", "getfasta", "-name", "-fi", reference, "-bed", ref_tes_gff, "-fo", ref_te_fasta]
        mccutils.run_command(command, log=log)
        mccutils.check_file_exists(ref_te_fasta)

    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...Failed to create TE fasta (bedtools getfasta) using reference:", reference, " and TE annotations:", ref_tes_gff, "check file formatting...exiting...", file=sys.stderr)
        sys.exit(1) 

    return ref_te_fasta

if __name__ == "__main__":
    main()