import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.jitterbug.jitterbug_post as config

def main():
    mccutils.log("jitterbug","jitterbug postprocessing")

    jitterbug_out = snakemake.input.jitterbug_out
    te_taxonomy = snakemake.input.taxonomy

    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")

    out = snakemake.output.out

    insertions = read_insertions(
                        jitterbug_out,
                        te_taxonomy,
                        chromosomes,
                        sample_name,
                        min_fwd_read_support=config.FILTER['MIN_FWD_READ_SUPPORT'],
                        min_rev_read_support=config.FILTER['MIN_REV_READ_SUPPORT'],
                        min_sr_support=config.FILTER['MIN_SPLIT_READ_SUPPORT'],
                        min_zygosity=config.FILTER['MIN_ZYGOSITY']
    )

    if len(insertions) >= 1:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="jitterbug")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="jitterbug")
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_nonredundant.bed"])

    # mccutils.run_command(["touch", out])


def read_insertions(jitterbug_gff, taxonomy, chroms, sample_name, min_fwd_read_support=0, min_rev_read_support=0, min_sr_support=0, min_zygosity=0.0):
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
                    sr = False
                    family = "NONE"
                    for feat in feats:
                        if "softclipped_pos" in feat:
                            pos = feat.split("=")[1]
                            pos = pos.replace("(","")
                            pos = pos.replace(")","")
                            pos = pos.split(",")
                            start = int(pos[0])-1
                            end = int(pos[1])

                            if start > -1 and end > -1:
                                insert.start = start
                                insert.end = end
                                sr = True
                        
                        if "predicted_superfam" in feat:
                            te  = feat.split("=")[1]
                            family = te_family[te]
                        
                        if "supporting_fwd_reads" in feat:
                            insert.jitterbug.supporting_fwd_reads = int(feat.split("=")[1])
                        
                        if "supporting_rev_reads" in feat:
                            insert.jitterbug.supporting_rev_reads = int(feat.split("=")[1])
                        
                        if "softclipped_support" in feat:
                            insert.jitterbug.split_read_support = int(feat.split("=")[1])
                        
                        if "zygosity" in feat:
                            insert.jitterbug.zygosity = float(feat.split("=")[1])
                
                    insert.name = family+"_non-reference_"+sample_name+"_jitterbug_"
                    if sr:
                        insert.name += "sr_"
                    else:
                        insert.name = "rp_"

                    if (
                        (insert.jitterbug.supporting_fwd_reads >= min_fwd_read_support) and 
                        (insert.jitterbug.supporting_rev_reads >= min_rev_read_support) and
                        (insert.jitterbug.split_read_support >= min_sr_support) and
                        (insert.jitterbug.zygosity >= min_zygosity)
                    ):
                        insertions.append(insert)
    
    return insertions



if __name__ == "__main__":                
    main()