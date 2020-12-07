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
    consensus = snakemake.input.consensus
    mcc_out = snakemake.params.mcc_out
    run_id = snakemake.params.run_id
    out_consensus = snakemake.output.consensus

    mccutils.log("processing","making consensus fasta")

    if not os.path.exists(mcc_out+"/tmp"):
        mccutils.mkdir(mcc_out+"/tmp")

    tmp = mcc_out+"/tmp/"+str(run_id)+"consensus.tmp"
    consensus = fix_fasta_lines(consensus, tmp)
    mccutils.replace_special_chars_fasta(consensus, out_consensus)
    

    mccutils.log("processing","consensus fasta created")


def fix_fasta_lines(infasta, outfasta, length=80):
    lines = fix_fasta.fix_fasta_lines(infasta, length)
    with open(outfasta, "w") as fa:
        for line in lines:
            fa.write(line+"\n")
    
    return outfasta





if __name__ == "__main__":                
    main()