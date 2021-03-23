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
    mccutils.log("popoolationte","processing PopoolationTE results")
    popoolationte_out = snakemake.input.popoolationte_out
    genome_fasta = snakemake.input.ref

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    chromosomes = snakemake.params.chromosomes.split(",")
    status_log = snakemake.params.status_log
    
    succeeded = mccutils.check_status_file(status_log)
    if succeeded:
        insertions = read_insertions(popoolationte_out, sample_name, chromosomes, require_both_end_support=config.REQUIRE_BOTH_END_SUPPORT, percent_read_support_threshold=config.PERCENT_READ_SUPPORT_THRESHOLD)
        if len(insertions) >= 1:
            insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="popoolationte")
            insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="popoolationte")
            output.write_vcf(insertions, genome_fasta, sample_name, "popoolationte", out_dir)
        else:
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_redundant.bed"])
            mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_nonredundant.bed"])
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_popoolationte_nonredundant.bed"])
    mccutils.log("popoolationte","PopoolationTE postprocessing complete")


def read_insertions(popoolationte, sample_name, chromosomes, require_both_end_support=True, percent_read_support_threshold=0.1):
    insertions = []

    with open(popoolationte, "r") as tsv:
        for line in tsv:
            insert = output.Insertion(output.Popoolationte())
            split_line = line.split("\t")
            insert.chromosome = split_line[0]
            pos_in_reference_seq = to_number(split_line[1])
            insert.support_info.support["flanks_supported"].value = split_line[2]
            insert.family = split_line[3]
            insert.support_info.support["frequency"].value = to_number(split_line[4], to_float=True)
            ref_te_id = split_line[6]
            insert.support_info.support["forward_insert_start"].value = to_number(split_line[8])
            insert.support_info.support["forward_insert_end"].value = to_number(split_line[9])
            insert.support_info.support["forward_insert_freq"].value = to_number(split_line[10], to_float=True)
            insert.support_info.support["forward_insert_cov"].value = to_number(split_line[11])
            insert.support_info.support["forward_presence_reads"].value = to_number(split_line[12])
            insert.support_info.support["forward_absence_reads"].value = to_number(split_line[13])
            insert.support_info.support["reverse_insert_start"].value = to_number(split_line[15])
            insert.support_info.support["reverse_insert_end"].value = to_number(split_line[16])
            insert.support_info.support["reverse_insert_freq"].value = to_number(split_line[17], to_float=True)
            insert.support_info.support["reverse_insert_cov"].value = to_number(split_line[18])
            insert.support_info.support["reverse_presence_reads"].value = to_number(split_line[19])
            insert.support_info.support["reverse_absence_reads"].value = to_number(split_line[20])
            
            if insert.support_info.support["forward_insert_start"].value == 0:
                insert.start = pos_in_reference_seq
                insert.end = insert.support_info.support["reverse_insert_start"].value
            
            elif insert.support_info.support["reverse_insert_start"].value == 0:
                insert.start = insert.support_info.support["forward_insert_end"].value
                insert.end = pos_in_reference_seq
            
            else:
                insert.start = insert.support_info.support["forward_insert_end"].value
                insert.end = insert.support_info.support["reverse_insert_start"].value           

            if "-" == ref_te_id:
                insert.type = "non-reference"
                insert.name = insert.family+"|non-reference|"+str(insert.support_info.support["frequency"].value)+"|"+sample_name+"|popoolationte|rp|"
            else:
                insert.type = "reference"
                insert.name = insert.family+"|reference|"+str(insert.support_info.support["frequency"].value)+"|"+sample_name+"|popoolationte|rp|"

            if not require_both_end_support:
                if ("FR" in insert.support_info.support["flanks_supported"].value and 
                        (insert.support_info.support["frequency"].value >= percent_read_support_threshold) and 
                        insert.chromosome in chromosomes):
                    insertions.append(insert)

                elif ("F" in insert.support_info.support["flanks_supported"].value and
                        (insert.support_info.support["forward_insert_freq"].value >= percent_read_support_threshold) and 
                        insert.chromosome in chromosomes):
                    insertions.append(insert)

                elif ((insert.support_info.support["reverse_insert_freq"].value >= percent_read_support_threshold) and
                        insert.chromosome in chromosomes):
                    insertions.append(insert)
            else:
                if ("FR" in insert.support_info.support["flanks_supported"].value and 
                        (insert.support_info.support["forward_insert_freq"].value >= percent_read_support_threshold or
                            insert.support_info.support["reverse_insert_freq"].value >= percent_read_support_threshold) and 
                        insert.chromosome in chromosomes):
                    insertions.append(insert)

    return insertions


def to_number(value, to_float=False):
    if not to_float:
        try:
            out_val = int(float(value))
        except ValueError as e:
            out_val = 0
    else:
        try:
            out_val = float(value)
        except ValueError as e:
            out_val = 0.0
    
    return out_val


if __name__ == "__main__":                
    main()