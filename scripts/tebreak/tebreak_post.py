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
    mccutils.log("tebreak","running tebreak post processing")
    tebreak_out = snakemake.input.tebreak_out
    ref_fasta = snakemake.input.ref_fasta

    out_dir = snakemake.params.out_dir
    ref_name = snakemake.params.ref_name
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    status_log = snakemake.params.status_log
    vcf_options = snakemake.params.vcf.split(",")

    prev_steps_succeeded = mccutils.check_status_file(status_log)
    if prev_steps_succeeded:
        insertions = read_insertions(tebreak_out, sample_name, chromosomes, config)

        if len(insertions) > 0:
            insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="tebreak")
            insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="tebreak")
            output.write_vcf(insertions, ref_fasta, sample_name, "tebreak", out_dir, vcf_options)
        else:
            mccutils.run_command(["touch", out_dir+"/"+sample_name+"_tebreak_redundant.bed"])
            mccutils.run_command(["touch", out_dir+"/"+sample_name+"_tebreak_nonredundant.bed"])
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_tebreak_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_tebreak_nonredundant.bed"])
    
    mccutils.log("tebreak","tebreak postprocessing complete")


def read_insertions(tebreak_out, sample_name, chromosomes, config):
    insertions = []
    header = {}
    with open(tebreak_out, "r") as inf:
        for ln,line in enumerate(inf):
            line = line.replace("\n","")
            split_line = line.split("\t")
            if ln == 0:
                for x,val in enumerate(split_line):
                    header[val] = x
            else:
                insert = output.Insertion(output.Tebreak())
                insert.chromosome = split_line[header['Chromosome']]
                insert.start = int(split_line[header['3_Prime_End']])+1
                insert.end = int(split_line[header['5_Prime_End']])
                insert.family = split_line[header['Superfamily']]
                insert.type = "non-reference"
                if split_line[header['Orient_5p']] == split_line[header['Orient_3p']]:
                    insert.strand = split_line[header['Orient_5p']]
                else:
                    insert.strand = "."
                
                if insert.strand == "-":
                    tmp = insert.start
                    insert.start = insert.end
                    insert.end = tmp

                insert.support_info.support["five_p_elt_match"].value = float(split_line[header['5p_Elt_Match']])
                insert.support_info.support["three_p_elt_match"].value = float(split_line[header['3p_Elt_Match']])
                insert.support_info.support["five_p_genome_match"].value = float(split_line[header['5p_Genome_Match']])
                insert.support_info.support["three_p_genome_match"].value = float(split_line[header['3p_Genome_Match']])
                insert.support_info.support["split_reads_5prime"].value = int(split_line[header['Split_reads_5prime']])
                insert.support_info.support["split_reads_3prime"].value = int(split_line[header['Split_reads_3prime']])
                insert.support_info.support["remapped_discordant"].value = int(split_line[header['Remapped_Discordant']])
                insert.support_info.support["remap_disc_fraction"].value = float(split_line[header['Remap_Disc_Fraction']])
                insert.support_info.support["remapped_splitreads"].value = int(split_line[header['Remapped_Splitreads']])
                insert.support_info.support["remap_split_fraction"].value = float(split_line[header['Remap_Split_Fraction']])

                insert.name = insert.family+"|non-reference|NA|"+sample_name+"|tebreak|sr|"

                if (
                    insert.chromosome in chromosomes and
                    insert.support_info.support["five_p_elt_match"].value >= config.MIN_5P_ELT_MATCH and 
                    insert.support_info.support["three_p_elt_match"].value >= config.MIN_3P_ELT_MATCH and
                    insert.support_info.support["five_p_genome_match"].value >= config.MIN_5P_GENOME_MATCH and
                    insert.support_info.support["three_p_genome_match"].value >= config.MIN_3P_GENOME_MATCH and
                    insert.support_info.support["split_reads_5prime"].value >= config.MIN_SPLIT_READS_5P and
                    insert.support_info.support["split_reads_3prime"].value >= config.MIN_SPLIT_READS_3P and
                    insert.support_info.support["remapped_discordant"].value >= config.MIN_REMAPPED_DISCORDANT and
                    insert.support_info.support["remap_disc_fraction"].value >= config.MIN_REMAP_DISC_FRACTION and
                    insert.support_info.support["remapped_splitreads"].value >= config.MIN_REMAPPED_SPLITREADS and
                    insert.support_info.support["remap_split_fraction"].value >= config.MIN_REMAP_SPLIT_FRACTION
                ):
                    insertions.append(insert)
    
    return insertions


if __name__ == "__main__":                
    main()