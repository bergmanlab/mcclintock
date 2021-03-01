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
    mccutils.log("retroseq","processing RetroSeq results")
    retroseq_out = snakemake.input.retroseq_out
    reference_fasta = snakemake.input.reference_fasta

    out_dir = snakemake.params.out_dir
    ref_name = snakemake.params.ref_name
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")

    insertions = read_insertions(retroseq_out, sample_name, chromosomes, support_threshold=config.READ_SUPPORT_THRESHOLD, breakpoint_threshold=config.BREAKPOINT_CONFIDENCE_THRESHOLD)
    if len(insertions) >= 1:
        insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="retroseq")
        insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="retroseq")
        output.write_vcf(insertions, reference_fasta, sample_name, "retroseq", out_dir)
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_retroseq_redundant.bed"])
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_retroseq_nonredundant.bed"])
    mccutils.log("retroseq","RetroSeq post processing complete")

def read_insertions(retroseq_vcf, sample_name, chromosomes, support_threshold=0, breakpoint_threshold=6):
    insertions = []

    with open(retroseq_vcf, "r") as vcf:
        for line in vcf:
            if "#" not in line:
                insert = output.Insertion(output.Retroseq())
                line = line.replace("\n","")
                split_line = line.split("\t")
                insert.chromosome = split_line[0]

                info = {}
                split_info = split_line[7].split(";")
                for i in split_info:
                    if "=" in i:
                        info[i.split("=")[0]] = i.split("=")[1]
                
                insert.family = (info['MEINFO'].split(",")[0]).split("-")[0]
                insert.start = int(info['MEINFO'].split(",")[1])
                insert.end = int(info['MEINFO'].split(",")[2])
                
                format_keys = split_line[8].split(":")
                format_vals = split_line[9].split(":")
                form = {}
                for x,key in enumerate(format_keys):
                    form[key] = format_vals[x]
                
                insert.support_info.support['spanning_pairs'].value = int(form['SP'])
                insert.support_info.support['supporting_reads'].value = int(form['GQ'])
                insert.support_info.support['clip3'].value = int(form['CLIP3'])
                insert.support_info.support['clip5'].value = int(form['CLIP5'])
                insert.support_info.support['call_status'].value = int(form['FL'])
                insert.type = "non-reference"
                insert.name = insert.family+"|non-reference|NA|"+sample_name+"|retroseq|rp|"

                if insert.support_info.support['supporting_reads'].value >= support_threshold and insert.support_info.support['call_status'].value >= breakpoint_threshold and insert.chromosome in chromosomes:
                    insertions.append(insert)
    
    return insertions


if __name__ == "__main__":                
    main()