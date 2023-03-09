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
    mccutils.log("temp","running TEMP post processing")
    insert_summary = snakemake.input.insert_summary
    absence_summary = snakemake.input.absence_summary
    te_gff = snakemake.input.te_gff
    reference_fasta = snakemake.input.reference_fasta
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    out_dir = snakemake.params.out_dir
    status_log = snakemake.params.status_log
    vcf_options = snakemake.params.vcf.split(",")

    prev_steps_succeeded = mccutils.check_status_file(status_log)

    if prev_steps_succeeded:
        insertions = read_insertion_summary(insert_summary, sample_name)
        absence_bed = make_absence_bed(absence_summary, sample_name, out_dir)
        non_absent_ref_insertions = get_non_absent_ref_tes(te_gff, absence_bed, sample_name, out_dir, log)
        insertions += non_absent_ref_insertions
        insertions = filter_insertions(insertions, chromosomes, acceptable_classes=config.PARAMS["acceptable_insertion_support_classes"], frequency_theshold=config.PARAMS["frequency_threshold"])
        if len(insertions) > 0:
            insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="temp")
            insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="temp")
            output.write_vcf(insertions, reference_fasta, sample_name, "temp", out_dir, vcf_options)
        else:
            mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp_redundant.bed"])
            mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp_nonredundant.bed"])
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp_nonredundant.bed"])
    mccutils.log("temp","TEMP postprocessing complete")






def read_insertion_summary(infile, sample):
    insertions = []
    with open(infile,"r") as inf:
        for x,line in enumerate(inf):
            if x > 0:
                insert = output.Insertion(output.Temp())
                split_line = line.split("\t")
                if len(split_line) == 14:
                    insert.chromosome = split_line[0]
                    insert.start = int(split_line[1])-1
                    insert.end = int(split_line[2])
                    insert.family = split_line[3]
                    insert.name = insert.family+"|non-reference|"+split_line[7]+"|"+sample+"|temp|"

                    if  "antisense" in split_line[4]:
                        insert.strand = "-"
                    else:
                        insert.strand = "+"
                        
                    insert.support_info.support['class'].value = split_line[5]
                    insert.support_info.support['variantsupport'].value = int(float(split_line[6]))
                    insert.support_info.support['frequency'].value = float(split_line[7])
                    insert.support_info.support['junction1'].value = int(split_line[8])
                    insert.support_info.support['junction1support'].value = int(split_line[9])
                    insert.support_info.support['junction2'].value = int(split_line[10])
                    insert.support_info.support['junction2support'].value = int(split_line[11])
                    insert.support_info.support['fiveprimesupport'].value = int(float(split_line[12]))
                    insert.support_info.support['threeprimesupport'].value = int(float(split_line[13].replace("\n","")))
                    insert.type = "non-reference"

                    if insert.end >= insert.start and insert.end > 0 and insert.start > -1:

                        # if split read, use junction positions as start and end
                        if insert.support_info.support['junction1support'].value > 0 and insert.support_info.support['junction2support'].value > 0:
                            insert.start = insert.support_info.support['junction1'].value
                            insert.end = insert.support_info.support['junction2'].value
                            insert.name = insert.name+"sr|"

                        # read pair
                        else:
                            insert.name = insert.name+"rp|" 

                        insertions.append(insert)
                    else:
                        print("<TEMP POST> Omitting malformed line from insertion summary results:", line)
                else:
                    print("<TEMP POST> Omitting malformed line from insertion summary results:", line)
    
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
                insert.name = split_line[9].split("=")[1]+"|reference|NA|"+sample+"|temp|nonab|"
                insert.strand = split_line[6]
                insert.type = "reference"
                
                insertions.append(insert)
    
    mccutils.remove(tmp_gff)

    return insertions
    

def filter_insertions(insertions, chromosomes, acceptable_classes=["1p1"], frequency_theshold=0.1):
    out = []
    for insert in insertions:
        if ( insert.chromosome in chromosomes and 
           ( insert.type == "reference" or 
               (insert.support_info.support['class'].value in acceptable_classes and 
                insert.support_info.support['frequency'].value > frequency_theshold))):
            out.append(insert)
    
    return out


if __name__ == "__main__":                
    main()