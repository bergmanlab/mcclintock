import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.popoolationte.popoolationte_post as config

class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "None"
        self.strand = "."
        self.support_type = ""
        self.f_read_support = 0
        self.r_read_support = 0
        self.f_read_support_percent = 0
        self.r_read_support_percent = 0

def main():
    mccutils.log("popoolationte","processing PopoolationTE results")
    popoolationte_out = snakemake.input.popoolationte_out

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    chromosomes = snakemake.params.chromosomes.split(",")

    insertions = read_insertions(popoolationte_out, sample_name, chromosomes, require_both_end_support=config.REQUIRE_BOTH_END_SUPPORT, percent_read_support_threshold=config.PERCENT_READ_SUPPORT_THRESHOLD)
    if len(insertions) >= 1:
        insertions = make_redundant_bed(insertions, sample_name, out_dir)
        make_nonredundant_bed(insertions, sample_name, out_dir)
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_nonredundant.bed"])
    mccutils.log("popoolationte","PopoolationTE postprocessing complete")


def read_insertions(popoolationte, sample_name, chromosomes, require_both_end_support=True, percent_read_support_threshold=0.1):
    insertions = []

    with open(popoolationte, "r") as tsv:
        for line in tsv:
            insert = Insertion()
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

            insert.support_type = split_line[2]

            if "F" in insert.support_type:
                insert.f_read_support = int(split_line[12])
                insert.f_read_support_percent = float(split_line[10])
            
            if "R" in insert.support_type:
                insert.r_read_support = int(split_line[19])
                insert.r_read_support_percent = float(split_line[17])

            if not require_both_end_support:
                if (insert.f_read_support_percent >= percent_read_support_threshold or insert.r_read_support_percent >= percent_read_support_threshold) and insert.chromosome in chromosomes:
                    insertions.append(insert)
            
            else:
                if "FR" in insert.support_type and (insert.f_read_support_percent >= percent_read_support_threshold or insert.r_read_support_percent >= percent_read_support_threshold) and insert.chromosome in chromosomes:
                    insertions.append(insert)
            
    return insertions

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

    redundant_bed = out_dir+"/"+sample_name+"_popoolationte_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_PoPoolationTE" description="'+sample_name+'_PoPoolationTE"\n'
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
            if (uniq_inserts[key].f_read_support + uniq_inserts[key].r_read_support) >  (insert.f_read_support + insert.r_read_support):
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

    nonredundant_bed = out_dir+"/"+sample_name+"_popoolationte_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_PoPoolationTE" description="'+sample_name+'_PoPoolationTE"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)

    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)



if __name__ == "__main__":                
    main()