import os
import sys
import subprocess
import importlib.util as il
spec = il.spec_from_file_location("config", snakemake.params.config)
config = il.module_from_spec(spec)
sys.modules[spec.name] = config
spec.loader.exec_module(config)
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import scripts.output as output



def main():
    relocate_gff = snakemake.input.relocate_gff
    te_gff = snakemake.input.te_gff
    reference_fasta = snakemake.input.reference_fasta

    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    status_log = snakemake.params.status_log
    vcf_options = snakemake.params.vcf.split(",")

    mccutils.log("relocate","processing RelocaTE results")

    prev_step_succeeded = mccutils.check_status_file(status_log)

    if prev_step_succeeded:
        insertions = get_insertions(
                        relocate_gff, 
                        sample_name, 
                        chromosomes, 
                        ref_l_threshold=config.PARAMS["ref_left_threshold"], 
                        ref_r_threshold=config.PARAMS["ref_right_threshold"], 
                        nonref_l_threshold=config.PARAMS["nonref_left_threshold"], 
                        nonref_r_threshold=config.PARAMS["nonref_right_threshold"]
                    )
        insertions = set_ref_orientations(insertions, te_gff)

        if len(insertions) >= 1:
            insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="relocate")
            insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="relocate")
            output.write_vcf(insertions, reference_fasta, sample_name, "relocate", out_dir, vcf_options)
        else:
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_relocate_redundant.bed"])
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_relocate_nonredundant.bed"])

    else:
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_relocate_redundant.bed"])
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_relocate_nonredundant.bed"])
    mccutils.log("relocate","RelocaTE postprocessing complete")


def get_insertions(gff, sample_name, chromosomes, ref_l_threshold=0, ref_r_threshold=0, nonref_l_threshold=0, nonref_r_threshold=0):
    insertions = []
    with open(gff, "r") as ingff:
        for line in ingff:
            if "#" not in line:
                split_line = line.split("\t")
                feats = split_line[8].split(";")
                insert = output.Insertion(output.Relocate())
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])
                insert.end = int(split_line[4])
                insert.strand = split_line[6]

                feat_id = ""
                feat_te_name = ""
                for feat in feats:
                    if "ID=" in feat:
                        feat_id = feat.split("=")[1]
                    elif "TE_Name=" in feat:
                        feat_te_name = feat.split("=")[1]
                    elif "Note=" in feat:
                        if "Shared" in feat:
                            insert.type = "reference"
                        elif "Non-reference" in feat:
                            insert.type = "non-reference"
                        else:
                            insert.type = "missing"
                    
                    elif "left_flanking_read_count=" in feat:
                        insert.support_info.support['left_flanking_reads'].value = int(feat.split("=")[1])
                    
                    elif "right_flanking_read_count=" in feat:
                        insert.support_info.support['right_flanking_reads'].value = int(feat.split("=")[1])
                
                if insert.type == "reference":
                    insert.family = feat_te_name
                    insert.name = feat_te_name+"|reference|NA|"+sample_name+"|relocate|sr|"
                elif insert.type == "non-reference":
                    feat_te_name = feat_id.split(".")[0]
                    insert.family = feat_te_name
                    insert.name = feat_te_name+"|non-reference|NA|"+sample_name+"|relocate|sr|"
            
            if insert.type == "reference" and insert.support_info.support['left_flanking_reads'].value >= ref_l_threshold and insert.support_info.support['right_flanking_reads'].value >= ref_r_threshold and insert.chromosome in chromosomes:
                insertions.append(insert)
            elif insert.type == "non-reference" and insert.support_info.support['left_flanking_reads'].value >= nonref_l_threshold and insert.support_info.support['right_flanking_reads'].value >= nonref_r_threshold and insert.chromosome in chromosomes:
                insertions.append(insert)
    
    return insertions

def set_ref_orientations(insertions, te_gff):
    ref_strands = {}
    out_inserts = []
    with open(te_gff, "r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                chrom = split_line[0]
                start = split_line[3]
                end = split_line[4]
                strand = split_line[6]
                ref_strands["_".join([chrom, start, end])] = strand
    
    for insert in insertions:
        if insert.type == "reference":
            insert.strand = ref_strands["_".join([insert.chromosome, str(insert.start), str(insert.end)])]
        
        out_inserts.append(insert)
    
    return out_inserts


if __name__ == "__main__":                
    main()