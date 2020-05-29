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
    masked_ref = snakemake.output.masked_ref


    ref_fasta = mask_reference(ref_fasta, te_gff, run_id, log, mcc_out)

    mccutils.run_command(["mv", ref_fasta, masked_ref])


# masks reference genome using reference TEs
def mask_reference(reference, ref_tes_gff, run_id, log, out):
    try:
        masked_reference = out+"/tmp/"+run_id+"tmpmaskedreference.fasta"
        command = ["bedtools", "maskfasta", "-fi", reference, "-fo", masked_reference, "-bed", ref_tes_gff]
        mccutils.run_command(command, log=log)

        mccutils.check_file_exists(masked_reference)
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...Failed to mask repeats (bedtools maskfasta) in reference: ", reference, " using repeat file:", ref_tes_gff, "check file formatting...exiting...", file=sys.stderr)
        sys.exit(1)

    return masked_reference

if __name__ == "__main__":
    main()