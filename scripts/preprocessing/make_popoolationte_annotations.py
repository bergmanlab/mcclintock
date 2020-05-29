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
    te_gff = snakemake.input.te_gff
    taxonomy = snakemake.input.taxonomy
    consensus = snakemake.input.consensus
    mcc_out = snakemake.params.mcc_out
    run_id = snakemake.params.run_id
    log = snakemake.params.log
    augment = snakemake.params.augment
    chromosomes = snakemake.params.chromosomes.split(",")
    popoolationte_taxonomy = snakemake.output.taxonomy
    popoolationte_te_gff = snakemake.output.te_gff
    
    if augment != "None":
        te_gff = remove_augmented_annotations(te_gff, chromosomes, mcc_out+"/tmp/"+str(run_id)+"ref_tes.gff")
        taxonomy = remove_augmented_taxonomy(taxonomy, te_gff, mcc_out+"/tmp/"+str(run_id)+"taxon.tsv")

    taxonomy = make_popoolationTE_taxonomy(taxonomy, consensus, run_id, mcc_out)

    mccutils.run_command(["cp", te_gff, popoolationte_te_gff])
    mccutils.run_command(["cp", taxonomy, popoolationte_taxonomy])




def remove_augmented_annotations(ingff, chromosomes, outgff):
    gff_lines = []
    with open(ingff,"r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                if split_line[0] in chromosomes:
                    gff_lines.append(line)
            else:
                gff_lines.append(line)
    
    with open(outgff,"w") as out:
        for line in gff_lines:
            out.write(line)
    
    return outgff

def remove_augmented_taxonomy(taxonomy, gff, outfile):
    tes_in_gff = []
    with open(gff, "r") as g:
        for line in g:
            if "#" not in line:
                split_line = line.split("\t")
                feats = split_line[8]
                split_feats = feats.split(";")
                te = split_feats[0]
                te = te.replace("ID=","")
                tes_in_gff.append(te)
    
    with open(taxonomy, "r") as taxon:
        with open(outfile,"w") as out:
            for line in taxon:
                split_line = line.split("\t")
                if split_line[0] in tes_in_gff:
                    out.write(line)
    
    return taxonomy

def make_popoolationTE_taxonomy(taxonomy, consensus, run_id, out):
    te_families = []
    popoolationTE_taxonomy = out+"/tmp/"+run_id+"tmppopoolationtetaxonomy.tsv"

    with open(taxonomy, "r") as intsv:
        for line in intsv:
            split_line = line.split("\t")
            te_families.append([split_line[0], split_line[1].replace("\n","")])
    
    fasta_records = SeqIO.parse(consensus,"fasta")

    for record in fasta_records:
        element = str(record.id)
        te_families.append([element, element])
    
    with open(popoolationTE_taxonomy, "w") as outtsv:
        header = "\t".join(["insert","id","family","superfamily","suborder","order","class","problem\n"])
        outtsv.write(header)

        for te_fam in te_families:
            element = te_fam[0]
            family = te_fam[1]
            line = "\t".join([element, family, family, family, "na", "na", "na", "0\n"])
            outtsv.write(line)
    
    return popoolationTE_taxonomy

if __name__ == "__main__":
    main()