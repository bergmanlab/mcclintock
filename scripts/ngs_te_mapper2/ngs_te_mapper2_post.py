import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import scripts.output as output
import config.ngs_te_mapper2.ngs_te_mapper2_post as config


def main():
    ref_bed = snakemake.input.ref_bed
    nonref_bed = snakemake.input.nonref_bed
    reference_fasta = snakemake.input.reference_fasta

    threads = snakemake.threads
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    out_dir = snakemake.params.out_dir
    chromosomes = snakemake.params.chromosomes.split(",")

    out_bed = snakemake.output[0]

    

    mccutils.log("ngs_te_mapper2","processing ngs_te_mapper2 results", log=log)

    insertions = read_insertions(ref_bed, nonref_bed, chromosomes, sample_name, out_dir)

    if len(insertions) > 0:
        insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="ngs_te_mapper2")
        intertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="ngs_te_mapper2")
        output.write_vcf(insertions, reference_fasta, sample_name, "ngs_te_mapper2", out_dir)

    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_ngs_te_mapper2_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_ngs_te_mapper2_nonredundant.bed"])

    mccutils.log("ngs_te_mapper","ngs_te_mapper postprocessing complete")

def read_insertions(ref_bed, nonref_bed, chromosomes, sample_name, out_dir):
    insertions = []
    with open(ref_bed,"r") as inbed:
        for line in inbed:
            line = line.replace("\n","")
            insert = output.Insertion(output.Ngs_te_mapper2())
            split_line = line.split("\t")
            insert.chromosome = split_line[0]
            insert.start = int(split_line[1])+1
            insert.end = int(split_line[2])
            insert.type = "reference"
            insert.strand = split_line[5]
            insert.family = split_line[3]
            insert.name = insert.family+"|"+insert.type+"|NA|"+sample_name+"|ngs_te_mapper2|sr|"
            if insert.chromosome in chromosomes:
                insertions.append(insert)
    
    with open(nonref_bed,"r") as inbed:
        for line in inbed:
            line = line.replace("\n","")
            insert = output.Insertion(output.Ngs_te_mapper2())
            split_line = line.split("\t")
            insert.chromosome = split_line[0]
            insert.start = int(split_line[1])+1
            insert.end = int(split_line[2])
            insert.type = "non-reference"
            insert.strand = split_line[5]
            insert.family = split_line[3].split("|")[0]
            insert.support_info.support['frequency'].value = float(split_line[3].split("|")[2])
            insert.support_info.support['three_prime_support'].value = int(split_line[3].split("|")[3])
            insert.support_info.support['five_prime_support'].value = int(split_line[3].split("|")[4])
            insert.support_info.support['reference_reads'].value = int(split_line[3].split("|")[5])
            insert.name = insert.family+"|"+insert.type+"|"+str(insert.support_info.support['frequency'].value)+"|"+sample_name+"|ngs_te_mapper2|sr|"
            if insert.chromosome in chromosomes:
                insertions.append(insert)

    return insertions


if __name__ == "__main__":                
    main()