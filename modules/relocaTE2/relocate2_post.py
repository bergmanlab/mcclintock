import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils
import config.relocate2.relocate2_post as config


class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "None"
        self.left_support = -1
        self.right_support = -1
        self.left_junction = -1
        self.right_junction = -1
        self.strand = "."

def main():
    nonref_gff = snakemake.input.nonref_gff
    ref_gff = snakemake.input.ref_gff
    rm_out = snakemake.input.rm_out

    log = snakemake.params.log
    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    

    ref_insertions = get_insertions(ref_gff, 
                                    sample_name, 
                                    insert_type="ref", 
                                    l_support_threshold=config.REF_LEFT_SUPPORT_THRESHOLD, 
                                    r_support_threshold=config.REF_RIGHT_SUPPORT_THRESHOLD,
                                    l_junction_threshold=config.REF_LEFT_JUNCTION_THRESHOLD,
                                    r_junction_threshold=config.REF_RIGHT_JUNCTION_THRESHOLD
                                    )

    nonref_insertions = get_insertions(nonref_gff, 
                                    sample_name, 
                                    insert_type="nonref", 
                                    l_support_threshold=config.NONREF_LEFT_SUPPORT_THRESHOLD, 
                                    r_support_threshold=config.NONREF_RIGHT_SUPPORT_THRESHOLD,
                                    l_junction_threshold=config.NONREF_LEFT_JUNCTION_THRESHOLD,
                                    r_junction_threshold=config.NONREF_RIGHT_JUNCTION_THRESHOLD
                                    )

    ref_insertions = fix_ref_te_names(ref_insertions, rm_out, sample_name)

    all_insertions = ref_insertions + nonref_insertions

    all_insertions = make_redundant_bed(all_insertions, sample_name, out_dir)

    make_nonredundant_bed(all_insertions, sample_name, out_dir)


def get_insertions(gff, sample_name, l_support_threshold=0, r_support_threshold=0, l_junction_threshold=0, r_junction_threshold=0, insert_type="ref"):
    insertions = []
    with open(gff, "r") as ingff:
        for line in ingff:
            if "#" not in line:
                line = line.replace(";","\t")
                split_line = line.split("\t")
                insert = Insertion()                  
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])
                insert.end = int(split_line[4])
                insert.strand = split_line[6]
                insert.type = insert_type

                insert.name = split_line[8].split("=")[1]

                if insert_type == "ref":
                    insert.right_junction = int(split_line[11].split(":")[1])
                    insert.left_junction = int(split_line[12].split(":")[1])
                    insert.right_support = int(split_line[13].split(":")[1])
                    insert.left_support = int(split_line[14].split(":")[1])
                else:
                    te_name = split_line[9].split("=")[1]
                    insert.name = te_name+"_non-reference_"+sample_name+"_relocate2_sr_"
                    insert.right_junction = int(split_line[12].split("=")[1])
                    insert.left_junction = int(split_line[13].split("=")[1])
                    insert.right_support = int(split_line[14].split("=")[1])
                    insert.left_support = int(split_line[15].split("=")[1])

                if insert.right_junction >= r_junction_threshold and insert.left_junction >= l_junction_threshold and insert.right_support >= r_support_threshold and insert.left_support >= l_support_threshold:
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
        insert.name = te_names[insert.name]+"_reference_"+sample_name+"_relocate2_sr_"
        out_insertions.append(insert)
    
    return out_insertions

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

    redundant_bed = out_dir+"/"+sample_name+"_relocate2_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_RelocaTE2" description="'+sample_name+'_RelocaTE2"\n'
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

    nonredundant_bed = out_dir+"/"+sample_name+"_relocate2_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_RelocaTE2" description="'+sample_name+'_RelocaTE2"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)
    

    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)

if __name__ == "__main__":                
    main()