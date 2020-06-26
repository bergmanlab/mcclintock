import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.telocate.telocate_post as config

class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "None"
        self.strand = "."
        self.read_pair_support = -1


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
        insertions = make_redundant_bed(insertions, sample_name, out_dir)
        make_nonredundant_bed(insertions, sample_name, out_dir)
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_telocate_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_telocate_nonredundant.bed"])
    mccutils.log("te-locate", "TE-Locate post processing complete")



def read_insertions(telocate_out, sample_name, chromosomes, rp_threshold=0):
    insertions = []
    with open(telocate_out,"r") as raw:
        for x, line in enumerate(raw):
            if x > 1:
                insert = Insertion()
                split_line = line.split("\t")
                insert.chromosome = split_line[0]
                insert.start = int(split_line[1])
                
                te_name = split_line[3].split("/")[1]
                if "old" in split_line[15]:
                    insert.type = "ref"
                    insert.end = insert.start+int(split_line[2])
                    insert.name = te_name+"_reference_"+sample_name+"_telocate_rp_"
                else:
                    insert.type = "nonref"
                    insert.end = insert.start
                    insert.name = te_name+"_non-reference_"+sample_name+"_telocate_rp_"

                if split_line[12] == "parallel":
                    insert.strand = "+"
                elif split_line[12] == "uncertain":
                    insert.strand = "."
                else:
                    insert.strand = "-"

                insert.read_pair_support = int(split_line[7])

                if insert.read_pair_support >= rp_threshold and insert.chromosome in chromosomes:
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
        if insert.type == "nonref":
            passed_insertions.append(insert)
        else:
            if insert.chromosome+"_"+str(insert.start) in gff_insertions.keys():
                insert.strand = gff_insertions[insert.chromosome+"_"+str(insert.start)]
                passed_insertions.append(insert)
            
    
    return passed_insertions


def make_redundant_bed(insertions, sample_name, out_dir):
    tmp_bed = out_dir+"/tmp_telocate_redundant.bed"

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

    redundant_bed = out_dir+"/"+sample_name+"_telocate_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_TE-locate" description="'+sample_name+'_TE-locate"\n'
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
            if uniq_inserts[key].read_pair_support >  insert.read_pair_support:
                uniq_inserts[key] = insert
    
    tmp_bed = out_dir+"/tmp_telocate_nonredundant.bed"
    with open(tmp_bed, "w") as outbed:
        for key in uniq_inserts.keys():
            insert = uniq_inserts[key]
            out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
            outbed.write(out_line+"\n")
    
    sorted_bed = out_dir+"/sorted.bed"
    command = ["bedtools", "sort", "-i", tmp_bed]
    mccutils.run_command_stdout(command, sorted_bed)

    nonredundant_bed = out_dir+"/"+sample_name+"_telocate_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_TE-locate" description="'+sample_name+'_TE-locate"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)
    

    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)


if __name__ == "__main__":                
    main()