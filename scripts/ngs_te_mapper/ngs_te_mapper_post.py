import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import scripts.output as output
import config.ngs_te_mapper.ngs_te_mapper_post as config


def main():
    raw_bed = snakemake.input.raw_bed
    reference_fasta = snakemake.input.reference_fasta

    threads = snakemake.threads
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    out_dir = snakemake.params.out_dir
    chromosomes = snakemake.params.chromosomes.split(",")

    out_bed = snakemake.output[0]

    

    mccutils.log("ngs_te_mapper","processing ngs_te_mapper results", log=log)

    insertions = read_insertions(raw_bed, chromosomes, sample_name, out_dir, min_read_cutoff=config.MIN_READ_SUPPORT)

    if len(insertions) > 0:
        insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="ngs_te_mapper")
        intertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="ngs_te_mapper")
        output.write_vcf(insertions, reference_fasta, sample_name, "ngs_te_mapper", out_dir)

    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_ngs_te_mapper_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_ngs_te_mapper_nonredundant.bed"])

    mccutils.log("ngs_te_mapper","ngs_te_mapper postprocessing complete")

def read_insertions(bed, chromosomes, sample_name, out_dir, min_read_cutoff=0):
    insertions = []
    with open(bed,"r") as inbed:
        for line in inbed:
            insert = output.Insertion(output.Ngs_te_mapper())
            line = line.replace(";","\t")
            split_line = line.split("\t")
            insert.chromosome = split_line[0]
            insert.start = int(split_line[1])+1
            insert.end = int(split_line[2])
            insert.type = split_line[8].replace("\n","")
            insert.strand = split_line[4]
            insert.family = split_line[5]
            insert.name = insert.family+"|"+insert.type+"|NA|"+sample_name+"|ngs_te_mapper|sr|"
            insert.support_info.support['supportingreads'].value = int(split_line[7])
            if insert.support_info.support['supportingreads'].value > min_read_cutoff and insert.chromosome in chromosomes:
                insertions.append(insert)

    return insertions


if __name__ == "__main__":                
    main()