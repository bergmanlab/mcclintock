import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.teflon.teflon_post as config

def main():
    mccutils.log("teflon","TEFLoN postprocessing")

    teflon_raw = snakemake.input.teflon_out
    ref_te_bed = snakemake.input.ref_bed
    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes
    out = snakemake.output.out

    ref_tes = get_ref_tes(ref_te_bed)
    insertions = read_insertions(
        teflon_raw, 
        chromosomes, 
        sample_name, 
        ref_tes,
        min_presence=config.PARAMETERS['min_presence_reads'], 
        max_absence=config.PARAMETERS['max_absence_reads'],
        min_presence_fraction=config.PARAMETERS['min_presence_fraction']
    )
    if len(insertions) >= 1:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="teflon")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="teflon")
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_teflon_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_teflon_nonredundant.bed"])

    mccutils.run_command(["touch", out])


def get_ref_tes(ref_te_bed):
    ref_tes = {}
    with open(ref_te_bed,"r") as bed:
        for line in bed:
            split_line = line.split("\t")
            ref_tes[split_line[3]] = [split_line[0], int(split_line[1]), int(split_line[2])]
    
    return ref_tes

def read_insertions(predictions, chroms, sample, ref_tes, min_presence=3, max_absence=None, min_presence_fraction=0.7):
    insertions = []

    with open(predictions, "r") as tsv:
        for line in tsv:
            split_line = line.split("\t")
            insert = mccutils.Insertion()

            insert.chromosome = split_line[0]

            if insert.chromosome in chroms:
                if split_line[1] != "-":
                    insert.start = int(split_line[1])
                else:
                    insert.start = int(split_line[2])-1
                
                if split_line[2] != "-":
                    insert.end = int(split_line[2])
                else:
                    insert.start = int(split_line[1])+1

                insert.family = split_line[3]
                
                insert.strand = split_line[5]
                
                # if reference prediction, uses ref TE coordinates
                if split_line[6] != "-":
                    insert.type = "reference"
                    insert.chromosome = ref_tes[split_line[6]][0]
                    insert.start = ref_tes[split_line[6]][1]
                    insert.end = ref_tes[split_line[6]][2]

                else:
                    insert.type = "non-reference"
                
                if split_line[7] == "+":
                    insert.teflon.left_sc_support = True
                
                if split_line[8] == "+":
                    insert.teflon.right_sc_support = True
                
                insert.teflon.presence_reads = int(split_line[9])
                insert.teflon.absence_reads = int(split_line[10])
                insert.teflon.ambiguous_reads = int(split_line[11])
                insert.teflon.allele_frequency = float(split_line[12])

                insert.name = split_line[3]+"_"+insert.type+"_"+sample+"_teflon_"

                if (
                    (insert.teflon.presence_reads >= min_presence) 
                    and (max_absence is None or insert.teflon.absence_reads <= max_absence)
                    and (insert.teflon.allele_frequency >= min_presence_fraction)
                ):
                    insertions.append(insert)
    
    return insertions
            



if __name__ == "__main__":                
    main()