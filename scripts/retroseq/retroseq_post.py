import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.retroseq.retroseq_post as config


def main():
    mccutils.log("retroseq","processing RetroSeq results")
    retroseq_out = snakemake.input.retroseq_out

    out_dir = snakemake.params.out_dir
    ref_name = snakemake.params.ref_name
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")

    insertions = read_insertions(retroseq_out, sample_name, chromosomes, support_threshold=config.READ_SUPPORT_THRESHOLD, breakpoint_threshold=config.BREAKPOINT_CONFIDENCE_THRESHOLD)
    if len(insertions) >= 1:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="retroseq")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="retroseq")
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_retroseq_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_retroseq_nonredundant.bed"])
    mccutils.log("retroseq","RetroSeq post processing complete")

def read_insertions(retroseq_vcf, sample_name, chromosomes, support_threshold=0, breakpoint_threshold=6):
    insertions = []

    with open(retroseq_vcf, "r") as vcf:
        for line in vcf:
            if "#" not in line:
                insert = mccutils.Insertion()
                line = line.replace(":","\t")
                line = line.replace("=", "\t")
                line = line.replace(",", "\t")
                split_line = line.split("\t")
                insert.chromosome = split_line[0]
                insert.start = int(split_line[10])
                insert.end = int(split_line[11])
                insert.retroseq.read_support = int(split_line[6])
                insert.name = split_line[9]+"_non-reference_"+sample_name+"_retroseq_rp_"
                insert.retroseq.breakpoint_confidence = int(split_line[22])

                if insert.retroseq.read_support >= support_threshold and insert.retroseq.breakpoint_confidence >= breakpoint_threshold and insert.chromosome in chromosomes:
                    insertions.append(insert)
    
    return insertions


if __name__ == "__main__":                
    main()