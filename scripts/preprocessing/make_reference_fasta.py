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
    reference = snakemake.input.ref
    augment = snakemake.params.augment
    mcc_out = snakemake.params.mcc_out
    run_id = snakemake.params.run_id
    log = snakemake.params.log
    out_ref = snakemake.output.ref
    out_aug_ref = snakemake.output.aug_ref

    if not os.path.exists(mcc_out+"/tmp"):
        mccutils.mkdir(mcc_out+"/tmp")

    tmp = mcc_out+"/tmp/"+str(run_id)+"reference.tmp"
    reference = fix_fasta_lines(reference, tmp)
    reference = mccutils.replace_special_chars_fasta(reference, tmp+"1")
    augmented_reference = reference
    if augment != "None":
        augment = fix_fasta_lines(augment, tmp+"2")
        augment = mccutils.replace_special_chars_fasta(augment, tmp+"3")
        augmented_reference = augment_reference(reference, augment, tmp+"4")
    
    mccutils.run_command(["cp", reference, out_ref])
    mccutils.run_command(["cp", augmented_reference, out_aug_ref])





def augment_reference(fasta1, fasta2, outfasta):
    lines = []
    with open(fasta1,"r") as fa1:
        for line in fa1:
            lines.append(line)
    
    with open(fasta2,"r") as fa2:
        for line in fa2:
            lines.append(line)
    
    with open(outfasta, "w") as out:
        for line in lines:
            out.write(line)
    
    return outfasta


def fix_fasta_lines(infasta, outfasta, length=80):
    lines = fix_fasta.fix_fasta_lines(infasta, length)
    with open(outfasta, "w") as fa:
        for line in lines:
            fa.write(line+"\n")
    
    return outfasta





if __name__ == "__main__":                
    main()