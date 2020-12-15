import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import scripts.output as output
import config.relocate2.relocate2_post as config


def main():
    nonref_gff = snakemake.input.nonref_gff
    ref_gff = snakemake.input.ref_gff
    rm_out = snakemake.input.rm_out

    log = snakemake.params.log
    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    
    mccutils.log("relocate2", "processing RelocaTE2 results")

    ref_insertions = get_insertions(ref_gff, 
                                    sample_name,
                                    chromosomes, 
                                    insert_type="ref", 
                                    l_support_threshold=config.REF_LEFT_SUPPORT_THRESHOLD, 
                                    r_support_threshold=config.REF_RIGHT_SUPPORT_THRESHOLD,
                                    l_junction_threshold=config.REF_LEFT_JUNCTION_THRESHOLD,
                                    r_junction_threshold=config.REF_RIGHT_JUNCTION_THRESHOLD)

    nonref_insertions = get_insertions(nonref_gff, 
                                    sample_name,
                                    chromosomes, 
                                    insert_type="nonref", 
                                    l_support_threshold=config.NONREF_LEFT_SUPPORT_THRESHOLD, 
                                    r_support_threshold=config.NONREF_RIGHT_SUPPORT_THRESHOLD,
                                    l_junction_threshold=config.NONREF_LEFT_JUNCTION_THRESHOLD,
                                    r_junction_threshold=config.NONREF_RIGHT_JUNCTION_THRESHOLD)

    ref_insertions = fix_ref_te_names(ref_insertions, rm_out, sample_name)

    all_insertions = ref_insertions + nonref_insertions

    if len(all_insertions) >= 1:
        all_insertions = output.make_redundant_bed(all_insertions, sample_name, out_dir, method="relocate2")
        output.make_nonredundant_bed(all_insertions, sample_name, out_dir, method="relocate2")
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_relocate2_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_relocate2_nonredundant.bed"])

    mccutils.log("relocate2", "RelocaTE2 postprocessing complete")

def get_insertions(gff, sample_name, chromosomes, l_support_threshold=0, r_support_threshold=0, l_junction_threshold=0, r_junction_threshold=0, insert_type="ref"):
    insertions = []
    with open(gff, "r") as ingff:
        for line in ingff:
            if "#" not in line:
                line = line.replace(";","\t")
                split_line = line.split("\t")
                insert = output.Insertion(output.Relocate2())                  
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])
                insert.end = int(split_line[4])
                insert.strand = split_line[6]
                insert.type = insert_type

                insert.name = split_line[8].split("=")[1]

                te_name = ""
                if insert_type == "ref":
                    insert.type = "reference"
                    insert.support_info.support['right_junction_reads'].value = int(split_line[11].split(":")[1])
                    insert.support_info.support['left_junction_reads'].value = int(split_line[12].split(":")[1])
                    insert.support_info.support['right_support_reads'].value = int(split_line[13].split(":")[1])
                    insert.support_info.support['left_support_reads'].value = int(split_line[14].split(":")[1])
                else:
                    insert.type = "non-reference"
                    te_name = split_line[9].split("=")[1]
                    te_name = te_name.split("/")[0]
                    insert.name = te_name+"|non-reference|NA|"+sample_name+"|relocate2|sr|"
                    insert.support_info.support['right_junction_reads'].value = int(split_line[12].split("=")[1])
                    insert.support_info.support['left_junction_reads'].value = int(split_line[13].split("=")[1])
                    insert.support_info.support['right_support_reads'].value = int(split_line[14].split("=")[1])
                    insert.support_info.support['left_support_reads'].value = int(split_line[15].split("=")[1])

                if ( insert.support_info.support['right_junction_reads'].value >= r_junction_threshold and 
                        insert.support_info.support['left_junction_reads'].value >= l_junction_threshold and 
                        insert.support_info.support['right_support_reads'].value >= r_support_threshold and 
                        insert.support_info.support['left_support_reads'].value >= l_support_threshold and 
                        insert.chromosome in chromosomes and
                        te_name != "repeat_name"):
                    insertions.append(insert)
    
    return insertions


def fix_ref_te_names(insertions, repeatmaskerout, sample_name):
    te_names = {}

    out_insertions = []
    with open(repeatmaskerout, "r") as infile:
        for x,line in enumerate(infile):
            if x > 2:
                tmp = line.split(" ")
                split_line = []
                for val in tmp:
                    if val != "":
                        split_line.append(val)
                
                chrom = split_line[4]
                start = split_line[5]
                end = split_line[6]
                te_name = split_line[9]

                key = "repeat_"+chrom+"_"+start+"_"+end
                te_names[key] = te_name
    
    for insert in insertions:
        insert.name = te_names[insert.name]+"|reference|NA|"+sample_name+"|relocate2|sr|"
        out_insertions.append(insert)
    
    return out_insertions


if __name__ == "__main__":                
    main()