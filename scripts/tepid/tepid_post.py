import os
import sys
import subprocess
import traceback
import importlib.util as il
spec = il.spec_from_file_location("config", snakemake.params.config)
config = il.module_from_spec(spec)
sys.modules[spec.name] = config
spec.loader.exec_module(config)
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import scripts.output as output

def main():
    insertions_bed = snakemake.input.insertions_bed
    deletions_bed = snakemake.input.deletions_bed
    insertions_support = snakemake.input.insertions_support
    deletions_support = snakemake.input.deletions_support
    te_gff = snakemake.input.te_gff
    te_taxonomy = snakemake.input.te_taxonomy
    reference_fasta = snakemake.input.reference_fasta

    chromosomes = snakemake.params.chromosomes.split(",")
    status_log = snakemake.params.status_log

    sample_name = snakemake.params.sample_name
    out_dir = snakemake.params.out_dir

    
    prev_steps_succeeded = mccutils.check_status_file(status_log)

    mccutils.log("tepid","running TEPID post processing")

    if prev_steps_succeeded:
        te_to_family = get_te_family_map(te_taxonomy)
        te_pos_to_family = get_te_pos_family_map(te_gff, te_to_family)
        insertions = read_insertions(insertions_bed, te_to_family, sample_name, te_pos_to_family, chromosomes, reference=False)
        insertions = add_support(insertions, insertions_support, threshold=config.READ_SUPPORT_THRESHOLD)

        deletions = read_insertions(deletions_bed, te_to_family, sample_name, te_pos_to_family, chromosomes, reference=True)
        deletions = add_support(deletions, deletions_support, threshold=config.READ_SUPPORT_THRESHOLD)
        non_abs_ref_insertions = get_non_absent_ref_tes(deletions, te_gff, te_to_family, sample_name)

        insertions += non_abs_ref_insertions
        if len(insertions) > 0:
            insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="tepid")
            insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="tepid")
            output.write_vcf(insertions, reference_fasta, sample_name, "tepid", out_dir)
        else:
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_tepid_redundant.bed"])
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_tepid_nonredundant.bed"])
            
    
    else:
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_tepid_redundant.bed"])
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_tepid_nonredundant.bed"])

    mccutils.log("tepid","TEPID post processing complete")
    

def get_te_family_map(taxonomy):
    te_to_family = {}
    with open(taxonomy,"r") as tsv:
        for line in tsv:
            split_line = line.split("\t")
            te_to_family[split_line[0]] = split_line[1].replace("\n","")
    
    return te_to_family

def get_te_pos_family_map(gff, te_to_family):
    pos_to_te = {}
    with open(gff, "r") as g:
        for line in g:
            if "#" not in line:
                split_line = line.split("\t")
                chrom = split_line[0]
                start = split_line[3]
                end = split_line[4]
                feats = split_line[8].split(";")
                te_name = ""
                for feat in feats:
                    if "ID=" in feat:
                        te_name = feat.split("=")[1].replace("\n","")
                pos_to_te[chrom+"_"+start+"_"+end] = te_to_family[te_name]
    
    return pos_to_te


def read_insertions(bed, te_to_family, sample_name, te_pos_to_family, chromosomes, reference=False):
    inserts = []
    with open(bed,"r") as b:
        for line in b:
            insert = output.Insertion(output.Tepid())
            split_line = line.split("\t")
            insert.chromosome = split_line[0]
            if insert.chromosome in chromosomes:
                insert.start = int(split_line[1])
                insert.end = int(split_line[2])

                if reference:
                    te_name = split_line[4].split(",")[0]
                    insert.family = te_to_family[te_name]
                    insert.strand = split_line[3]
                    insert.type = "reference"
                    insert.name = insert.family+"|reference|NA|"+sample_name+"|tepid|nonab|"
                else:
                    te_chrom = split_line[3]
                    te_start = split_line[4]
                    te_end = split_line[5]
                    insert.family = te_pos_to_family[te_chrom+"_"+te_start+"_"+te_end]
                    insert.type = "non-reference"
                    insert.name = insert.family+"|non-reference|NA|"+sample_name+"|tepid|sr|"
                
                insert.support_info.id = split_line[-1].replace("\n","")
                inserts.append(insert)
    
    return inserts
            

def add_support(inserts, support_file, threshold=0):
    filtered_inserts = []
    support = {}

    with open(support_file, "r") as txt:
        for line in txt:
            split_line = line.split("\t")
            support_id = split_line[0].replace(">","")
            support_val = len(split_line[1].split(","))
            support[support_id] = support_val
    
    for insert in inserts:
        insert.support_info.support['supporting_reads'].value = support[insert.support_info.id]
        if insert.support_info.support['supporting_reads'].value >= threshold:
            filtered_inserts.append(insert)
    

    return filtered_inserts


def get_non_absent_ref_tes(deletions, te_gff, te_to_family, sample_name):
    ref_tes = []
    with open(te_gff, "r") as gff:
        for line in gff:
            ref_te = output.Insertion(output.Tepid())
            split_line = line.split("\t")
            ref_te.chromosome = split_line[0]
            ref_te.start = int(split_line[3])
            ref_te.end = int(split_line[4])
            ref_te.strand = split_line[6]
            feats = split_line[8]
            split_feats = feats.split(";")
            te_id = ""
            for f in split_feats:
                if "ID=" in f:
                    te_id = f.split("=")[1]
            
            ref_te.family = te_to_family[te_id]
            ref_te.type = "reference"
            ref_te.name = ref_te.family+"|reference|NA|"+sample_name+"|tepid|nonab|"
            ref_tes.append(ref_te)
    
    absent = []
    for deletion in deletions:
        key = "_".join([deletion.chromosome, str(deletion.start), str(deletion.end), deletion.strand, deletion.family])
        absent.append(key)
    
    non_absent = []
    for te in ref_tes:
        key = "_".join([te.chromosome, str(te.start), str(te.end), te.strand, te.family])
        if key not in absent:
            non_absent.append(te)
    
    return non_absent




if __name__ == "__main__":                
    main()