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
    mccutils.log("jitterbug","jitterbug postprocessing")

    jitterbug_out = snakemake.input.jitterbug_out
    te_taxonomy = snakemake.input.taxonomy
    reference_fasta = snakemake.input.reference_fasta

    out_dir = snakemake.params.out_dir
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    status_log = snakemake.params.status_log
    vcf_options = snakemake.params.vcf.split(",")

    out = snakemake.output.out

    prev_steps_succeeded = mccutils.check_status_file(status_log)

    if prev_steps_succeeded:
        insertions = read_insertions(
                            jitterbug_out,
                            te_taxonomy,
                            chromosomes,
                            sample_name,
                            min_fwd_read_support=config.FILTER['MIN_FWD_READ_SUPPORT'],
                            min_rev_read_support=config.FILTER['MIN_REV_READ_SUPPORT'],
                            min_sr_support=config.FILTER['MIN_SPLIT_READ_SUPPORT'],
                            min_zygosity=config.FILTER['MIN_ZYGOSITY']
        )

        if len(insertions) >= 1:
            insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="jitterbug")
            insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="jitterbug")
            output.write_vcf(insertions, reference_fasta, sample_name, "jitterbug", out_dir, vcf_options)
        else:
            mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_redundant.bed"])
            mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_nonredundant.bed"])
    
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_jitterbug_nonredundant.bed"])

    # mccutils.run_command(["touch", out])


def read_insertions(jitterbug_gff, taxonomy, chroms, sample_name, min_fwd_read_support=0, min_rev_read_support=0, min_sr_support=0, min_zygosity=0.0):
    insertions = []

    te_family = {}
    with open(taxonomy,"r") as tsv:
        for line in tsv:
            line = line.replace("\n","")
            split_line = line.split("\t")
            te_family[split_line[0]] = split_line[1]

    with open(jitterbug_gff, "r") as gff:
        for line in gff:
            line = line.replace("\n","")
            split_line = line.split("\t")
            if len(split_line) == 9:
                # insert = mccutils.Insertion()
                insert = output.Insertion(output.Jitterbug())

                insert.chromosome = split_line[0]
                if insert.chromosome in chroms:
                    insert.start = int(split_line[3])
                    insert.end = int(split_line[4])
                    insert.type = "non-reference"

                    feats = split_line[8]
                    feats = feats.replace(" ","")
                    feats = feats.split(";")
                    supporting_families = []
                    sr = False
                    family = "NONE"
                    for feat in feats:
                        if "softclipped_pos" in feat:
                            pos = feat.split("=")[1]
                            pos = pos.replace("(","")
                            pos = pos.replace(")","")
                            pos = pos.split(",")
                            start = int(pos[0])-1
                            end = int(pos[1])

                            if start > -1 and end > -1:
                                insert.start = start
                                insert.end = end
                                sr = True
                        
                        if "predicted_superfam" in feat:
                            te  = feat.split("=")[1]
                            family = te_family[te]
                            insert.family = family
                        
                        if "supporting_fwd_reads" in feat:
                            insert.support_info.support['supporting_fwd_reads'].value = int(feat.split("=")[1])
                        
                        if "supporting_rev_reads" in feat:
                            insert.support_info.support['supporting_rev_reads'].value = int(feat.split("=")[1])
                        
                        if "softclipped_support" in feat:
                            insert.support_info.support['softclipped_support'].value = int(feat.split("=")[1])
                        
                        if "zygosity" in feat:
                            insert.support_info.support['zygosity'].value = float(feat.split("=")[1])
                
                    insert.name = family+"|non-reference|"+str(insert.support_info.support['zygosity'].value)+"|"+sample_name+"|jitterbug|"
                    if sr:
                        insert.name += "sr|"
                    else:
                        insert.name = "rp|"

                    if (
                        (insert.support_info.support['supporting_fwd_reads'].value >= min_fwd_read_support) and 
                        (insert.support_info.support['supporting_rev_reads'].value >= min_rev_read_support) and
                        (insert.support_info.support['softclipped_support'].value >= min_sr_support) and
                        (insert.support_info.support['zygosity'].value >= min_zygosity)
                    ):
                        insertions.append(insert)
    
    return insertions



if __name__ == "__main__":                
    main()