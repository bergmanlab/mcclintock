import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.popoolationte2.popoolationte2_post as config

class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "nonref"
        self.family = ""
        self.strand = "."
        self.support_type = ""
        self.frequency = -1
        self.added = False

def main():
    print("<POPOOLATIONTE2> Processing PopoolationTE2 results...")
    te_predictions = snakemake.input.popoolationte2_out
    te_gff = snakemake.input.te_gff
    taxonomy = snakemake.input.taxonomy
    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    log = snakemake.params.log

    ref_tes = get_ref_tes(te_gff, taxonomy, chromosomes)
    insertions = read_insertions(te_predictions, ref_tes, chromosomes, sample_name, both_end_support_needed=config.REQUIRE_BOTH_END_SUPPORT, support_threshold=config.FREQUENCY_THRESHOLD)
    insertions = make_redundant_bed(insertions, sample_name, out_dir)
    make_nonredundant_bed(insertions, sample_name, out_dir)

def get_ref_tes(gff, taxon, chroms):
    ref_inserts = []
    te_family = {}
    with open(taxon,"r") as t:
        for line in t:
            split_line = line.split("\t")
            te_id = split_line[0]
            family = split_line[1]
            te_family[te_id] = family

    with open(gff, "r") as g:
        for line in g:
            if "#" not in line:
                split_line = line.split("\t")
                insert = Insertion()
                insert.type = "ref"
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])
                insert.end = int(split_line[4])
                insert.strand = split_line[6]
                insert.family = te_family[split_line[2]]
                if insert.chromosome in chroms:
                    ref_inserts.append(insert)
    
    return ref_inserts

def read_insertions(predictions, ref_tes, chroms, sample, both_end_support_needed=True, support_threshold=0.1):
    insertions = []

    with open(predictions,"r") as tsv:
        for line in tsv:
            split_line = line.split("\t")
            insert = Insertion()
            insert.chromosome = split_line[1]
            insert.start = int(split_line[2])
            insert.end = int(split_line[2])
            insert.strand = split_line[3]
            insert.family = split_line[4]
            insert.support_type = split_line[6]
            insert.frequency = float(split_line[8])

            if (insert.support_type == "FR" or not both_end_support_needed) and insert.frequency > support_threshold:     
                # determine if insert is a ref insert       
                for x in range(0,len(ref_tes)):
                    if ref_tes[x].start <= insert.start and insert.start <= ref_tes[x].end:
                        insert.family = ref_tes[x].family
                        insert.added = ref_tes[x].added
                        if not ref_tes[x].added:
                            ref_tes[x].added = True
                            
                        insert.type = "ref"
                        insert.start = ref_tes[x].start
                        insert.end = ref_tes[x].end
                        insert.strand = ref_tes[x].strand
                
                if insert.type == "ref":
                    insert.name = insert.family+"_reference_"+str(insert.frequency)+"_"+sample+"_popoolationte2_rp_"
                else:
                    insert.name = insert.family+"_non-reference_"+str(insert.frequency)+"_"+sample+"_popoolationte2_rp_"
                
                if not insert.added:
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

    redundant_bed = out_dir+"/"+sample_name+"_popoolationte2_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_PoPoolationTE2" description="'+sample_name+'_PoPoolationTE2"\n'
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

    nonredundant_bed = out_dir+"/"+sample_name+"_popoolationte2_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_PoPoolationTE2" description="'+sample_name+'_PoPoolationTE2"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)

    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)


if __name__ == "__main__":                
    main()