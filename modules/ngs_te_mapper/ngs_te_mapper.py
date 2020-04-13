import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils
import config.ngs_te_mapper.ngs_te_mapper as config


def main():
    consensus_fasta = snakemake.input.consensus_fasta
    reference_fasta = snakemake.input.reference_fasta
    fastq1 = snakemake.input.fastq1
    fastq2 = snakemake.input.fastq2
    threads = snakemake.threads
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    script_dir = snakemake.params.script_dir
    out_dir = snakemake.params.out_dir

    is_paired = True
    if snakemake.params.raw_fq2 == "None":
        is_paired = False
    
    mccutils.mkdir(out_dir+"/unfiltered/")
    command = ['Rscript', "--vanilla", script_dir+"/ngs_te_mapper.R", "genome="+reference_fasta, "teFile="+consensus_fasta, "tsd="+str(config.MAX_TSD), "thread="+str(threads), "output="+out_dir+"/unfiltered/", "sourceCodeFolder="+script_dir]

    if is_paired:
        command.append("sample="+fastq1+";"+fastq2)
    else:
        command.append("sample="+fastq1)
    
    print("<NGS_TE_MAPPER> Running ngs_te_mapper...")
    mccutils.run_command(command, log=log)
    print("<NGS_TE_MAPPER> Processing ngs_te_mapper results...")

    fq1_basename = mccutils.get_base_name(fastq1)
    fq2_basename = mccutils.get_base_name(fastq2)
    if is_paired:
        raw_bed = out_dir+"/unfiltered/bed_tsd/"+fq1_basename+"_"+fq2_basename+"insertions.bed"
    else:
        raw_bed = out_dir+"/unfiltered/bed_tsd/"+fq1_basename+"insertions.bed"
    
    process_bed(raw_bed, sample_name, log, out_dir)

    print("<NGS_TE_MAPPER> ngs_te_mapper complete...")

    

def process_bed(bed, sample_name, log, out_dir):
    unsorted_bed = out_dir+"/unsorted.bed"
    with open(unsorted_bed, "w") as outbed:
        with open(bed,"r") as inbed:
            for x,line in enumerate(inbed):
                line = line.replace(";","\t")
                split_line = line.split("\t")
                outline = "\t".join([split_line[0], split_line[1], split_line[2], split_line[5]+"_"+split_line[8].replace("\n","")+"_"+sample_name+"_ngs_te_mapper_sr_"+str(x+1),"0", split_line[4]])
                outbed.write(outline+"\n")
    
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

    mccutils.remove(unsorted_bed)
    mccutils.remove(sorted_bed)

if __name__ == "__main__":                
    main()