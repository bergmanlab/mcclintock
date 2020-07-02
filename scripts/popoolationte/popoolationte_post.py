import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.popoolationte.popoolationte_post as config


def main():
    mccutils.log("popoolationte","processing PopoolationTE results")
    popoolationte_out = snakemake.input.popoolationte_out

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    chromosomes = snakemake.params.chromosomes.split(",")

    insertions = read_insertions(popoolationte_out, sample_name, chromosomes, require_both_end_support=config.REQUIRE_BOTH_END_SUPPORT, percent_read_support_threshold=config.PERCENT_READ_SUPPORT_THRESHOLD)
    if len(insertions) >= 1:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="popoolationte")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="popoolationte")
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_nonredundant.bed"])
    mccutils.log("popoolationte","PopoolationTE postprocessing complete")


def read_insertions(popoolationte, sample_name, chromosomes, require_both_end_support=True, percent_read_support_threshold=0.1):
    insertions = []

    with open(popoolationte, "r") as tsv:
        for line in tsv:
            insert = mccutils.Insertion()
            split_line = line.split("\t")
            insert.chromosome = split_line[0]
            if "-" in split_line[8]:
                insert.start = int(float(split_line[1]))
                insert.end = int(float(split_line[15]))
            
            elif "-" in split_line[15]:
                insert.start = int(float(split_line[9]))
                insert.end = int(float(split_line[1]))
            
            else:
                insert.start = int(float(split_line[9]))
                insert.end = int(float(split_line[15]))               

            if "-" in split_line[6]:
                insert.name = split_line[3]+"_non-reference_"+split_line[4]+"_"+sample_name+"_popoolationte_rp_"
            else:
                insert.name = split_line[3]+"_reference_"+split_line[4]+"_"+sample_name+"_popoolationte_rp_"

            insert.popoolationte.support_type = split_line[2]

            if "F" in insert.popoolationte.support_type:
                insert.popoolationte.f_read_support = int(split_line[12])
                insert.popoolationte.f_read_support_percent = float(split_line[10])
            
            if "R" in insert.popoolationte.support_type:
                insert.popoolationte.r_read_support = int(split_line[19])
                insert.popoolationte.r_read_support_percent = float(split_line[17])

            if not require_both_end_support:
                if (insert.popoolationte.f_read_support_percent >= percent_read_support_threshold or insert.popoolationte.r_read_support_percent >= percent_read_support_threshold) and insert.chromosome in chromosomes:
                    insertions.append(insert)
            
            else:
                if "FR" in insert.popoolationte.support_type and (insert.popoolationte.f_read_support_percent >= percent_read_support_threshold or insert.popoolationte.r_read_support_percent >= percent_read_support_threshold) and insert.chromosome in chromosomes:
                    insertions.append(insert)
            
    return insertions



if __name__ == "__main__":                
    main()