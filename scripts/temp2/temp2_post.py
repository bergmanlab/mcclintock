import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import scripts.output as output
import config.temp2.temp2_post as config

def main():
    mccutils.log("temp2","running TEMP2 post processing")
    insert_bed = snakemake.input.insert_bed
    absence_summary = snakemake.input.absence_summary
    te_gff = snakemake.input.te_gff
    reference_fasta = snakemake.input.reference_fasta
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    out_dir = snakemake.params.out_dir

    insertions = read_insertions(insert_bed, sample_name, chromosomes, config)
    absence_bed = make_absence_bed(absence_summary, sample_name, out_dir)
    non_absent_ref_insertions = get_non_absent_ref_tes(te_gff, absence_bed, sample_name, out_dir, log)
    insertions += non_absent_ref_insertions

    if len(insertions) > 0:
        insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="temp2")
        insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="temp2")
        output.write_vcf(insertions, reference_fasta, sample_name, "temp2", out_dir)
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp2_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp2_nonredundant.bed"])

    mccutils.log("temp2","TEMP2 postprocessing complete")


def read_insertions(insert_bed, sample_name, chromosomes, config):
    insertions = []
    with open(insert_bed, "r") as inf:
        for x,line in enumerate(inf):
            if x > 0:
                insert = output.Insertion(output.Temp2())
                split_line = line.split("\t")
                if len(split_line) == 15:
                    insert.chromosome = split_line[0]
                    insert.start = int(split_line[1])
                    insert.end = int(split_line[2])
                    insert.family = split_line[3].split(":")[0]
                    insert.type = "non-reference"
                    insert.support_info.support["frequency"].value = float(split_line[4])
                    insert.strand = split_line[5]
                    insert.support_info.support["class"].value = split_line[6]
                    insert.support_info.support["supportreads"].value = float(split_line[7])
                    insert.support_info.support["referencereads"].value = float(split_line[8])
                    insert.support_info.support["fiveprimesupport"].value = float(split_line[9])
                    insert.support_info.support["threeprimesupport"].value = float(split_line[10])
                    insert.support_info.support["reliability"].value = float(split_line[12].replace("%","")) # rare enties have a % sign for some reason
                    insert.support_info.support["fiveprimejunctionsupport"].value = float(split_line[13])
                    insert.support_info.support["threeprimejunctionsupport"].value = float(split_line[14])

                    insert.name = insert.family+"|non-reference|"+str(insert.support_info.support['frequency'].value)+"|"+sample_name+"|temp2|"

                    if (
                        insert.chromosome in chromosomes and
                        insert.support_info.support["frequency"].value >= config.FREQUENCY_THRESHOLD and 
                        insert.support_info.support["class"].value in config.ACCEPTABLE_INSERTION_SUPPORT_CLASSES
                    ):
                        insertions.append(insert)
    
    return insertions

def make_absence_bed(summary_file, sample, out):
    out_bed = out+"/"+sample+".absent.bed"
    lines = []
    with open(summary_file, "r") as inf:
        for x,line in enumerate(inf):
            if x > 0:
                split_line = line.split("\t")
                new_line = "\t".join([split_line[0], split_line[1], split_line[2]])
                new_line += "\n"
                lines.append(new_line)
    
    if len(lines) < 1:
        lines.append("empty\t0\t1\n")
    
    with open(out_bed,"w") as bed:
        for line in lines:
            bed.write(line)
    
    return out_bed

def get_non_absent_ref_tes(te_gff, absence_bed, sample, out, log):
    insertions = []
    tmp_gff = out+"/tmp.ref_nonabs.gff"
    command = ["bedtools", "subtract", "-A", "-a", te_gff, "-b", absence_bed]
    mccutils.run_command_stdout(command, tmp_gff, log=log)

    with open(tmp_gff,"r") as gff:
        for line in gff:
            if "#" not in line:
                line = line.replace(";","\t")
                split_line = line.split("\t")
                insert = output.Insertion(output.Temp())
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])
                insert.end = int(split_line[4])
                insert.name = split_line[9].split("=")[1]+"|reference|NA|"+sample+"|temp2|nonab|"
                insert.strand = split_line[6]
                insert.type = "reference"
                
                insertions.append(insert)
    
    mccutils.remove(tmp_gff)

    return insertions


if __name__ == "__main__":                
    main()