import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.relocate.relocate_post as config

class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "None"
        self.left_support = -1
        self.right_support = -1
        self.strand = "."

def main():
    relocate_gff = snakemake.input.relocate_gff
    te_gff = snakemake.input.te_gff

    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name

    print("<RELOCATE POST> Processing RelocaTE results...")

    insertions = get_insertions(relocate_gff, sample_name, ref_l_threshold=config.REF_LEFT_THRESHOLD, ref_r_threshold=config.REF_RIGHT_THRESHOLD, nonref_l_threshold=config.NONREF_LEFT_THRESHOLD, nonref_r_threshold=config.NONREF_RIGHT_THRESHOLD)

    insertions = set_ref_orientations(insertions, te_gff)

    insertions = make_redundant_bed(insertions, sample_name, out_dir)

    make_nonredundant_bed(insertions, sample_name, out_dir)

    print("<RELOCATE POST> RelocaTE postprocessing results...")


def get_insertions(gff, sample_name, ref_l_threshold=0, ref_r_threshold=0, nonref_l_threshold=0, nonref_r_threshold=0):
    insertions = []
    with open(gff, "r") as ingff:
        for line in ingff:
            if "#" not in line:
                split_line = line.split("\t")
                feats = split_line[8].split(";")
                insert = Insertion()
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
                            insert.type = "ref"
                        elif "Non-reference" in feat:
                            insert.type = "nonref"
                        else:
                            insert.type = "missing"
                    
                    elif "left_flanking_read_count=" in feat:
                        insert.left_support = int(feat.split("=")[1])
                    
                    elif "right_flanking_read_count=" in feat:
                        insert.right_support = int(feat.split("=")[1])
                
                if insert.type == "ref":
                    insert.name = feat_te_name+"_reference_"+sample_name+"_relocate_sr_"
                elif insert.type == "nonref":
                    feat_te_name = feat_id.split(".")[0]
                    insert.name = feat_te_name+"_non-reference_"+sample_name+"_relocate_sr_"
            
            if insert.type == "ref" and insert.left_support >= ref_l_threshold and insert.right_support >= ref_r_threshold:
                insertions.append(insert)
            elif insert.type == "nonref" and insert.left_support >= nonref_l_threshold and insert.right_support >= nonref_r_threshold:
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
        if insert.type == "ref":
            insert.strand = ref_strands["_".join([insert.chromosome, str(insert.start), str(insert.end)])]
        
        out_inserts.append(insert)
    
    return out_inserts


def make_redundant_bed(insertions, sample_name, out_dir):
    tmp_bed = out_dir+"/tmp.bed"

    insertion_dict = {}
    out_inserts = []
    for insert in insertions:
        insertion_dict[ "_".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])] = insert


    with open(tmp_bed, "w") as out:
        for insert in insertions:
            out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
            out.write(out_line+"\n")
    
    sorted_bed = out_dir+"/sorted.bed"
    command = ["bedtools", "sort", "-i", tmp_bed]
    mccutils.run_command_stdout(command, sorted_bed)

    redundant_bed = out_dir+"/"+sample_name+"_relocate_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_RelocaTE" description="'+sample_name+'_RelocaTE"\n'
        outbed.write(header)
        with open(sorted_bed, "r") as inbed:
            for x, line in enumerate(inbed):

                # outputs inserts in sorted order with unique number added to name
                key = line.replace("\t","_")
                key = key.replace("\n","")
                insert = insertion_dict[key]
                insert.name += str(x+1)
                out_inserts.append(insert)

                # write to bed with unique number added to name
                split_line = line.split("\t")
                split_line[3] += str(x+1)
                line = "\t".join(split_line)
                outbed.write(line)
    
    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)

    return out_inserts


def make_nonredundant_bed(insertions, sample_name, out_dir):
    uniq_inserts = {}

    for insert in insertions:
        key = "_".join([insert.chromosome, str(insert.end)])
        if key not in uniq_inserts.keys():
            uniq_inserts[key] = insert
        else:
            if (uniq_inserts[key].left_support + uniq_inserts[key].right_support) < (insert.left_support + insert.right_support):
                uniq_inserts[key] = insert
    
    tmp_bed = out_dir+"/tmp.bed"
    with open(tmp_bed, "w") as outbed:
        for key in uniq_inserts.keys():
            insert = uniq_inserts[key]
            out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
            outbed.write(out_line+"\n")
    
    sorted_bed = out_dir+"/sorted.bed"
    command = ["bedtools", "sort", "-i", tmp_bed]
    mccutils.run_command_stdout(command, sorted_bed)

    nonredundant_bed = out_dir+"/"+sample_name+"_relocate_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_RelocaTE" description="'+sample_name+'_RelocaTE"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)
    

    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)

if __name__ == "__main__":                
    main()