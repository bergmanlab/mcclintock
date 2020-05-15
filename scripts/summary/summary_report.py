import os
import sys
import subprocess
from Bio import SeqIO
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.config as config
from datetime import datetime


def main():
    out_files = snakemake.input.out_files
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2
    ref = snakemake.input.ref
    taxonomy = snakemake.input.taxonomy

    bam = snakemake.params.bam
    flagstat = snakemake.params.flagstat
    median_insert_size = snakemake.params.median_insert_size
    methods = snakemake.params.methods
    results_dir = snakemake.params.results_dir
    sample_name = snakemake.params.sample_name
    command = snakemake.params.command
    execution_dir = snakemake.params.execution_dir
    start_time = snakemake.params.time
    raw_fq2 = snakemake.params.raw_fq2
    out_dir = snakemake.params.out_dir

    paired = True
    if raw_fq2 == "None":
        paired = False

    out_file_map = {}
    for method in config.OUT_PATHS.keys():
        out_file_map[method] = config.OUT_PATHS[method]
        out_file_map[method] = out_file_map[method].replace("{{results}}", results_dir)
        out_file_map[method] = out_file_map[method].replace("{{samplename}}", sample_name)

    # make_te_summary_table(methods, out_file_map, snakemake.output.te_summary)
    make_te_csv(methods, out_file_map, taxonomy, snakemake.output.te_csv)
    make_run_summary(out_file_map, methods, fq1, fq2, ref, bam, flagstat, median_insert_size, command, execution_dir, start_time, out_dir, snakemake.output.summary_report, paired=paired)
    


def make_te_summary_table(methods, out_file_map, out_table):

    len_longest_name = 0
    for method in config.ALL_METHODS:
        if len(method) > len_longest_name:
            len_longest_name = len(method)


    width1 = len_longest_name+2
    width2 = 10
    width3 = 13
    width4 = 14

    with open(out_table,"w") as out:
        out.write("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")
        out.write(pad("METHOD", width1) + pad("ALL",width2) + pad("REFERENCE",width3) + pad("NON-REFERENCE", width4) + "\n")
        out.write("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")
        for method in config.ALL_METHODS:
            if "nonredundant.bed" in out_file_map[method]:
                if method in methods:
                    all_te, ref_te, nonref_te = get_te_counts(out_file_map[method])
                    out.write(pad(method, width1) + pad(str(all_te), width2) + pad(str(ref_te), width3) + pad(str(nonref_te), width4) + "\n")
                else:
                    out.write(pad(method, width1) + pad("NA", width2) + pad("NA", width3) + pad("NA", width4) + "\n")

        out.write("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")


def get_te_counts(bed):
    all_te = 0
    ref_te = 0
    nonref_te = 0

    with open(bed,"r") as b:
        for line in b:
            split_line = line.split("\t")
            if len(split_line) == 6:
                all_te += 1
                if "non-reference" in split_line[3]:
                    nonref_te += 1
                else:
                    ref_te += 1


    return all_te, ref_te, nonref_te


def pad(string, total_len, symbol=" ", front=False):
    if not front:
        string = string+(symbol*(total_len - len(string)))
    else:
        string = (symbol*(total_len - len(string)))+string

    return string
    


def make_te_csv(methods, out_file_map, taxonomy, out_csv):
    te_counts = {}
    method_counts = {}
    header = ["TE-Family"]

    # get all TE family names
    te_names = []
    with open(taxonomy,"r") as taxon:
        for line in taxon:
            te_name = line.split("\t")[1].replace("\n","")
            if te_name not in te_names:
                te_names.append(te_name)
    
    # setup TE family count dicts and header
    for x,te_name in enumerate(te_names):
        te_counts[te_name] = {}
        for method in config.ALL_METHODS:
            if "nonredundant.bed" in out_file_map[method]:
                if method in methods:
                    te_counts[te_name][method] = [0,0]
                else:
                    te_counts[te_name][method] = ["NA","NA"]
                
                if x < 1:
                    header += [method+"_all", method+"_ref", method+"_nonref"]
    
    
    # count TE family predictions for each method
    for method in methods:
        bed = out_file_map[method]
        if "nonredundant.bed" in bed:
            with open(bed, "r") as b:
                for line in b:
                    split_line = line.split("\t")
                    if len(split_line) == 6:
                        name = split_line[3]
                        te_name = ""
                        split_name = name.split("_")
                        for x,part in enumerate(split_name):
                            if part == "reference" or part == "non-reference":
                                te_name = "_".join(split_name[:x])
                        
                        if te_name == "":
                            print("ERROR: couldn't find TE name in:", name)
                    
                        else:
                            if te_name not in te_counts.keys():
                                te_counts[te_name] = {method:[0,0]}
                            elif method not in te_counts[te_name].keys():
                                te_counts[te_name][method] = [0,0]
                                
                            if "non-reference" in name:
                                te_counts[te_name][method][1] += 1
                            else:
                                te_counts[te_name][method][0] += 1

    # create te family count csv
    with open(out_csv,"w") as out:
        out.write(",".join(header)+"\n")
        for te in te_names:
            line = [te.replace(",","_")]
            for method in config.ALL_METHODS:
                if "nonredundant.bed" in out_file_map[method]:
                    if method in methods:
                        line += [str(te_counts[te][method][0]+te_counts[te][method][1]), str(te_counts[te][method][0]), str(te_counts[te][method][1])]
                    else:
                        line +=["NA","NA","NA"]
            
            out.write(",".join(line)+"\n")
    

def make_run_summary(out_file_map, methods, fq1, fq2, ref, bam, flagstat, median_insert_size, command, execution_dir, start_time, out_dir, out_file, paired=False):
    out_lines = []
    out_lines.append(("-"*34)+"\n")
    out_lines.append("MCCLINTOCK SUMMARY REPORT\n")
    out_lines.append(("-"*34) + "\n")
    split_command = command.split(" ")
    for x,split in enumerate(split_command):
        if split[0] == "-":
            split_command[x] = split_command[x].replace("-", " \\\n\t-",1)
    command = " ".join(split_command)
    out_lines.append("Command:\n"+command+"\n")
    out_lines.append("\nrun from directory: "+execution_dir+"\n")
    out_lines.append(pad("Started:",12)+start_time+"\n")
    out_lines.append(pad("Completed:",12)+"{{END_TIME}}"+"\n\n")

    if os.path.exists(bam) and os.path.exists(flagstat) and os.path.exists(median_insert_size):
        out_lines.append(("-"*34)+"\n")
        out_lines.append("MAPPED READ INFORMATION\n")
        out_lines.append(("-"*34) + "\n")

        out_lines.append(pad("read1 sequence length:",24) + str(estimate_read_length(fq1))+"\n")
        if paired:
            out_lines.append(pad("read2 sequence length:",24) + str(estimate_read_length(fq2))+"\n")
        
        with open(flagstat,"r") as stat:
            for line in stat:
                if "read1" in line:
                    reads = line.split("+")[0].replace(" ","")
                    out_lines.append(pad("read1 reads:",24) + str(reads) + "\n")
                
                if paired and "read2" in line:
                    reads = line.split("+")[0].replace(" ","")
                    out_lines.append(pad("read2 reads:",24) + str(reads) + "\n")

            
        with open(median_insert_size,"r") as median_insert:
            for line in median_insert:
                line = line.split("=")[1].replace("\n","")
                out_lines.append(pad("median insert size:",24) + str(int(float(line))) + "\n")
        
        out_lines.append(pad("avg genome coverage:",24) + str(get_avg_coverage(ref, bam, out_dir)) + "\n")
        out_lines.append(("-"*34)+"\n")


    len_longest_name = 0
    for method in config.ALL_METHODS:
        if len(method) > len_longest_name:
            len_longest_name = len(method)


    width1 = len_longest_name+2
    width2 = 10
    width3 = 13
    width4 = 14

    out_lines.append("\n")
    out_lines.append("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")
    out_lines.append(pad("METHOD", width1) + pad("ALL",width2) + pad("REFERENCE",width3) + pad("NON-REFERENCE", width4) + "\n")
    out_lines.append("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")
    for method in config.ALL_METHODS:
        if "nonredundant.bed" in out_file_map[method]:
            if method in methods:
                all_te, ref_te, nonref_te = get_te_counts(out_file_map[method])
                out_lines.append(pad(method, width1) + pad(str(all_te), width2) + pad(str(ref_te), width3) + pad(str(nonref_te), width4) + "\n")
            else:
                out_lines.append(pad(method, width1) + pad("NA", width2) + pad("NA", width3) + pad("NA", width4) + "\n")

    out_lines.append("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")

    
    with open(out_file,"w") as out:
        for line in out_lines:
            if "{{END_TIME}}" in line:
                now = datetime.now()
                completed = now.strftime("%Y-%m-%d %H:%M:%S")
                line = line.replace("{{END_TIME}}",completed)
            print(line, end="")
            out.write(line)
        



def estimate_read_length(fq, reads=100):
    lengths = []
    with open(fq,"r") as f:
        for x, line in enumerate(f):
            if x%4 == 1:
                lengths.append(len(line.replace('\n',"")))
            
            if x >= reads:
                break
    
    length = sum(lengths)//len(lengths)

    return length


def get_avg_coverage(ref, bam, out):
    chrom = []
    fasta_records = SeqIO.parse(ref,"fasta")
    for record in fasta_records:
        chrom.append(str(record.id))
    
    tmp = out+"/tmp"
    command = ['samtools','depth',bam]
    mccutils.run_command_stdout(command, tmp)

    cov_total = 0
    pos = 0 
    with open(tmp,"r") as depth:
        for line in depth:
            split_line = line.split("\t")
            if split_line[0] in chrom:
                pos += 1
                cov_total += int(split_line[2])
    
    mccutils.remove(tmp)
    return round(cov_total/pos, 3)

    


if __name__ == "__main__":                
    main()

