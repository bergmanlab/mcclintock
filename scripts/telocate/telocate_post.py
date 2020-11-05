import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.telocate.telocate_post as config


def main():
    mccutils.log("te-locate","processing TE-Locate results")
    telocate_raw = snakemake.input.telocate_raw
    te_gff = snakemake.input.te_gff

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")

    insertions = read_insertions(telocate_raw, sample_name, chromosomes, rp_threshold=config.READ_PAIR_SUPPORT_THRESHOLD)
    insertions = filter_by_reference(insertions, te_gff)
    if len(insertions) > 0:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="telocate")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir,method="telocate")
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_telocate_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_telocate_nonredundant.bed"])
    mccutils.log("te-locate", "TE-Locate post processing complete")



def read_insertions(telocate_out, sample_name, chromosomes, rp_threshold=0):
    insertions = []
    with open(telocate_out,"r") as raw:
        for x, line in enumerate(raw):
            if x > 1:
                insert = mccutils.Insertion()
                split_line = line.split("\t")
                insert.chromosome = split_line[0]
                insert.start = int(split_line[1])
                
                te_name = split_line[3].split("/")[1]
                if "old" in split_line[15]:
                    insert.type = "reference"
                    insert.end = insert.start+int(split_line[2])
                    insert.name = te_name+"|reference|"+sample_name+"|telocate|rp|"
                else:
                    insert.type = "non-reference"
                    insert.end = insert.start
                    insert.name = te_name+"|non-reference|"+sample_name+"|telocate|rp|"

                if split_line[12] == "parallel":
                    insert.strand = "+"
                elif split_line[12] == "uncertain":
                    insert.strand = "."
                else:
                    insert.strand = "-"

                insert.telocate.read_pair_support = int(split_line[7])

                if insert.telocate.read_pair_support >= rp_threshold and insert.chromosome in chromosomes:
                    insertions.append(insert)
    
    return insertions


def filter_by_reference(insertions, te_gff):
    gff_insertions = {}
    passed_insertions = []
    with open(te_gff,"r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                chrom = split_line[0]
                start = split_line[3]
                strand = split_line[6]
                gff_insertions[chrom+"_"+start] = strand

    
    for insert in insertions:
        if insert.type == "non-reference":
            passed_insertions.append(insert)
        else:
            if insert.chromosome+"_"+str(insert.start) in gff_insertions.keys():
                insert.strand = gff_insertions[insert.chromosome+"_"+str(insert.start)]
                passed_insertions.append(insert)
            
    
    return passed_insertions



if __name__ == "__main__":                
    main()