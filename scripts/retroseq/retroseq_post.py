import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.retroseq.retroseq_post as config

class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "None"
        self.strand = "."
        self.read_support = -1
        self.breakpoint_confidence = -1

def main():
    mccutils.log("retroseq","processing RetroSeq results")
    retroseq_out = snakemake.input.retroseq_out

    out_dir = snakemake.params.out_dir
    ref_name = snakemake.params.ref_name
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")

    insertions = read_insertions(retroseq_out, sample_name, chromosomes, support_threshold=config.READ_SUPPORT_THRESHOLD, breakpoint_threshold=config.BREAKPOINT_CONFIDENCE_THRESHOLD)
    if len(insertions) >= 1:
        insertions = make_redundant_bed(insertions, sample_name, out_dir)
        make_nonredundant_bed(insertions, sample_name, out_dir)
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_retroseq_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_retroseq_nonredundant.bed"])
    mccutils.log("retroseq","RetroSeq post processing complete")

def read_insertions(retroseq_vcf, sample_name, chromosomes, support_threshold=0, breakpoint_threshold=6):
    insertions = []

    with open(retroseq_vcf, "r") as vcf:
        for line in vcf:
            if "#" not in line:
                insert = Insertion()
                line = line.replace(":","\t")
                line = line.replace("=", "\t")
                line = line.replace(",", "\t")
                split_line = line.split("\t")
                insert.chromosome = split_line[0]
                insert.start = int(split_line[10])
                insert.end = int(split_line[11])
                insert.read_support = int(split_line[6])
                insert.name = split_line[9]+"_non-reference_"+sample_name+"_retroseq_rp_"
                insert.breakpoint_confidence = int(split_line[20])

                if insert.read_support >= support_threshold and insert.breakpoint_confidence >= breakpoint_threshold and insert.chromosome in chromosomes:
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

    redundant_bed = out_dir+"/"+sample_name+"_retroseq_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_RetroSeq" description="'+sample_name+'_RetroSeq"\n'
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
            if uniq_inserts[key].read_support >  insert.read_support:
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

    nonredundant_bed = out_dir+"/"+sample_name+"_retroseq_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_RetroSeq" description="'+sample_name+'_RetroSeq"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)

    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)


if __name__ == "__main__":                
    main()