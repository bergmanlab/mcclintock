import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
from collections import Counter

def main():
    mccutils.log("jitterbug","jitterbug postprocessing")

    jitterbug_out = snakemake.input.jitterbug_out
    te_taxonomy = snakemake.input.taxonomy

    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")

    out = snakemake.output.out

    insertions = read_insertions(jitterbug_out, te_taxonomy, chromosomes, sample_name)
    print("INSERTIONS:", len(insertions))

    if len(insertions) >= 1:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="jitterbug")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="jitterbug")
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_nonredundant.bed"])

    # mccutils.run_command(["touch", out])


def read_insertions(jitterbug_gff, taxonomy, chroms, sample_name):
    insertions = []

    te_family = {}
    with open(taxonomy,"r") as tsv:
        for line in tsv:
            line = line.replace("\n","")
            split_line = line.split("\t")
            te_family[split_line[0]] = split_line[1]

    with open(jitterbug_gff, "r") as gff:
        for line in gff:
            line = line.replace("\n","")
            split_line = line.split("\t")
            if len(split_line) == 9:
                insert = mccutils.Insertion()

                insert.chromosome = split_line[0]
                if insert.chromosome in chroms:
                    insert.start = int(split_line[3])
                    insert.end = int(split_line[4])
                    insert.type = "non-reference"

                    feats = split_line[8]
                    feats = feats.replace(" ","")
                    feats = feats.split(";")
                    supporting_families = []
                    for feat in feats:
                        if "softclipped_pos" in feat:
                            pos = feat.split("=")[1]
                            pos = pos.replace("(","")
                            pos = pos.replace(")","")
                            pos = pos.split(",")
                            start = int(pos[0])
                            end = int(pos[1])

                            if start > -1 and end > -1:
                                insert.start = start
                                insert.end = end
                        
                        if "Inserted_TE_tags" in feat:
                            te_list = feat.split("=")[1]
                            te_list = te_list.split(",")
                            for te in te_list:
                                family = te_family[te]
                                supporting_families.append(family)

                    family_counts = Counter(supporting_families)
                    best_support = ["None",0]
                    for fam in supporting_families:
                        if family_counts[fam] > best_support[1]:
                            best_support = [fam, family_counts[fam]]

                    insert.name = best_support[0]+"_non-reference_"+sample_name+"_jitterbug_"
                    insertions.append(insert)
    
    return insertions



if __name__ == "__main__":                
    main()