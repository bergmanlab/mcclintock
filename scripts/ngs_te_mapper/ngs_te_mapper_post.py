import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.ngs_te_mapper.ngs_te_mapper_post as config


def main():
    raw_bed = snakemake.input.raw_bed

    threads = snakemake.threads
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    out_dir = snakemake.params.out_dir
    chromosomes = snakemake.params.chromosomes.split(",")

    out_bed = snakemake.output[0]

    

    print("<NGS_TE_MAPPER POST> Processing ngs_te_mapper results...")

    process_bed(raw_bed, chromosomes, sample_name, log, out_dir, min_read_cutoff=config.MIN_READ_SUPPORT)

def process_bed(bed, chromosomes, sample_name, log, out_dir, min_read_cutoff=0):
    unsorted_bed = out_dir+"/unsorted.bed"
    with open(unsorted_bed, "w") as outbed:
        with open(bed,"r") as inbed:
            insertion_count = 0
            for x,line in enumerate(inbed):
                line = line.replace(";","\t")
                split_line = line.split("\t")
                if int(split_line[7]) > min_read_cutoff and split_line[0] in chromosomes:
                    insertion_count += 1
                    outline = "\t".join([split_line[0], split_line[1], split_line[2], split_line[5]+"_"+split_line[8].replace("\n","")+"_"+sample_name+"_ngs_te_mapper_sr_"+str(x+1),"0", split_line[4]])
                    outbed.write(outline+"\n")
    
    if insertion_count >= 1:
        sorted_bed = out_dir+"/sorted.bed"
        command = ["bedtools", "sort", "-i", unsorted_bed]
        mccutils.run_command_stdout(command, sorted_bed, log=log)

        final_bed = out_dir+"/"+sample_name+"_ngs_te_mapper_nonredundant.bed"
        with open(final_bed,"w") as outbed:
            header = 'track name="'+sample_name+'_ngs_te_mapper" description="'+sample_name+'_ngs_te_mapper"\n'
            outbed.write(header)
            with open(sorted_bed, "r") as inbed:
                for line in inbed:
                    # line = line.replace("NA",".")
                    outbed.write(line)
        mccutils.remove(sorted_bed)
    else:
        mccutils.run_command(["touch",out_dir+"/"+sample_name+"_ngs_te_mapper_nonredundant.bed"])
        
    mccutils.remove(unsorted_bed)



if __name__ == "__main__":                
    main()