import os
import sys
import subprocess
import traceback
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.tepid.tepid_post as config


def main():
    insertions_bed = snakemake.input.insertions_bed
    deletions_bed = snakemake.input.deletions_bed
    insertions_support = snakemake.input.insertions_support
    deletions_support = snakemake.input.deletions_support
    te_gff = snakemake.input.te_gff
    te_taxonomy = snakemake.input.te_taxonomy

    sample_name = snakemake.params.sample_name
    out_dir = snakemake.params.out_dir

    mccutils.log("tepid","running TEPID post processing")
    te_to_family = get_te_family_map(te_taxonomy)
    insertions = read_insertions(insertions_bed, te_to_family, sample_name, reference=False)
    insertions = add_support(insertions, insertions_support, threshold=config.READ_SUPPORT_THRESHOLD)

    deletions = read_insertions(deletions_bed, te_to_family, sample_name, reference=True)
    deletions = add_support(deletions, deletions_support, threshold=config.READ_SUPPORT_THRESHOLD)
    non_abs_ref_insertions = get_non_absent_ref_tes(deletions, te_gff, te_to_family, sample_name)

    insertions += non_abs_ref_insertions
    if len(insertions) > 0:
        mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="tepid")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="tepid")
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_tepid_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_tepid_nonredundant.bed"])
        
    mccutils.log("tepid","TEPID post processing complete")

def get_te_family_map(taxonomy):
    te_to_family = {}
    with open(taxonomy,"r") as tsv:
        for line in tsv:
            split_line = line.split("\t")
            te_to_family[split_line[0]] = split_line[1].replace("\n","")
    
    return te_to_family

def read_insertions(bed, te_to_family, sample_name, reference=False):
    inserts = []
    with open(bed,"r") as b:
        for line in b:
            insert = mccutils.Insertion()
            split_line = line.split("\t")
            insert.chromosome = split_line[0]
            insert.start = int(split_line[1])
            insert.end = int(split_line[2])

            if reference:
                te_name = split_line[4].split(",")[0]
                insert.family = te_to_family[te_name]
                insert.strand = split_line[3]
                insert.type = "reference"
                insert.name = insert.family+"_reference_"+sample_name+"_tepid_nonab_"
            else:
                te_name = split_line[6].split(",")[0]
                insert.family = te_to_family[te_name]
                insert.type = "non-reference"
                insert.name = insert.family+"_non-reference_"+sample_name+"_tepid_"
            
            insert.tepid.id = split_line[-1].replace("\n","")
            inserts.append(insert)
    
    return inserts
            

def add_support(inserts, support_file, threshold=0):
    filtered_inserts = []
    support = {}

    with open(support_file, "r") as txt:
        for line in txt:
            split_line = line.split("\t")
            support_id = split_line[0].replace(">","")
            support_val = len(split_line[1].split(","))
            support[support_id] = support_val
    
    for insert in inserts:
        insert.tepid.support = support[insert.tepid.id]
        if insert.tepid.support >= threshold:
            filtered_inserts.append(insert)
    

    return filtered_inserts


def get_non_absent_ref_tes(deletions, te_gff, te_to_family, sample_name):
    ref_tes = []
    with open(te_gff, "r") as gff:
        for line in gff:
            ref_te = mccutils.Insertion()
            split_line = line.split("\t")
            ref_te.chromosome = split_line[0]
            ref_te.start = int(split_line[3])
            ref_te.end = int(split_line[4])
            ref_te.strand = split_line[6]
            feats = split_line[8]
            split_feats = feats.split(";")
            te_id = ""
            for f in split_feats:
                if "ID=" in f:
                    te_id = f.split("=")[1]
            
            ref_te.family = te_to_family[te_id]
            ref_te.type = "reference"
            ref_te.name = ref_te.family+"_reference_"+sample_name+"_tepid_nonab_"
            ref_tes.append(ref_te)
    
    absent = []
    for deletion in deletions:
        key = "_".join([deletion.chromosome, str(deletion.start), str(deletion.end), deletion.strand, deletion.family])
        absent.append(key)
    
    non_absent = []
    for te in ref_tes:
        key = "_".join([te.chromosome, str(te.start), str(te.end), te.strand, te.family])
        if key not in absent:
            non_absent.append(te)
    
    return non_absent




if __name__ == "__main__":                
    main()