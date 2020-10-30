import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.relocate.relocate_post as config


def main():
    relocate_gff = snakemake.input.relocate_gff
    te_gff = snakemake.input.te_gff

    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")

    mccutils.log("relocate","processing RelocaTE results")

    insertions = get_insertions(relocate_gff, sample_name, chromosomes, ref_l_threshold=config.REF_LEFT_THRESHOLD, ref_r_threshold=config.REF_RIGHT_THRESHOLD, nonref_l_threshold=config.NONREF_LEFT_THRESHOLD, nonref_r_threshold=config.NONREF_RIGHT_THRESHOLD)

    insertions = set_ref_orientations(insertions, te_gff)

    if len(insertions) >= 1:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="relocate")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="relocate")
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_relocate_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_relocate_nonredundant.bed"])

    mccutils.log("relocate","RelocaTE postprocessing complete")


def get_insertions(gff, sample_name, chromosomes, ref_l_threshold=0, ref_r_threshold=0, nonref_l_threshold=0, nonref_r_threshold=0):
    insertions = []
    with open(gff, "r") as ingff:
        for line in ingff:
            if "#" not in line:
                split_line = line.split("\t")
                feats = split_line[8].split(";")
                insert = mccutils.Insertion()
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])
                insert.end = int(split_line[4])
                insert.strand = split_line[6]

                feat_id = ""
                feat_te_name = ""
                for feat in feats:
                    if "ID=" in feat:
                        feat_id = feat.split("=")[1]
                    elif "TE_Name=" in feat:
                        feat_te_name = feat.split("=")[1]
                    elif "Note=" in feat:
                        if "Shared" in feat:
                            insert.type = "reference"
                        elif "Non-reference" in feat:
                            insert.type = "non-reference"
                        else:
                            insert.type = "missing"
                    
                    elif "left_flanking_read_count=" in feat:
                        insert.relocate.left_support = int(feat.split("=")[1])
                    
                    elif "right_flanking_read_count=" in feat:
                        insert.relocate.right_support = int(feat.split("=")[1])
                
                if insert.type == "reference":
                    insert.name = feat_te_name+"|reference|"+sample_name+"|relocate|sr|"
                elif insert.type == "non-reference":
                    feat_te_name = feat_id.split(".")[0]
                    insert.name = feat_te_name+"|non-reference|"+sample_name+"|relocate|sr|"
            
            if insert.type == "reference" and insert.relocate.left_support >= ref_l_threshold and insert.relocate.right_support >= ref_r_threshold and insert.chromosome in chromosomes:
                insertions.append(insert)
            elif insert.type == "non-reference" and insert.relocate.left_support >= nonref_l_threshold and insert.relocate.right_support >= nonref_r_threshold and insert.chromosome in chromosomes:
                insertions.append(insert)
    
    return insertions

def set_ref_orientations(insertions, te_gff):
    ref_strands = {}
    out_inserts = []
    with open(te_gff, "r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                chrom = split_line[0]
                start = split_line[3]
                end = split_line[4]
                strand = split_line[6]
                ref_strands["_".join([chrom, start, end])] = strand
    
    for insert in insertions:
        if insert.type == "reference":
            insert.strand = ref_strands["_".join([insert.chromosome, str(insert.start), str(insert.end)])]
        
        out_inserts.append(insert)
    
    return out_inserts


if __name__ == "__main__":                
    main()