import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    consensus_fasta = snakemake.input.consensus_fasta
    te_gff = snakemake.input.te_gff
    reference_fasta = snakemake.input.reference_fasta
    fq1 = snakemake.input.fq1
    sq2 = snakemake.input.fq2

    log = snakemake.params.log
    raw_fq2 = snakemake.params.raw_fq2
    script_dir = snakemake.params.script_dir
    out_dir = snakemake.params.out_dir

    out_gff = snakemake.output[0]

    is_paired = True
    if snakemake.params.raw_fq2 == "None":
        is_paired = False

    annotation = make_annotation_file(te_gff, out_dir)
    os.chdir(out_dir)

    fq_dir = os.path.dirname(fq1)
    command = ["perl", script_dir+"/relocaTE.pl", "-t", consensus_fasta, "-d", fq_dir, "-g", reference_fasta, "-o", ".", "-r", annotation]

    if is_paired:
        command += ["-1", "_1", "-2", "_2"]
    else:
        command += ["-u", "unPaired"]
    
    print("<RELOCATE> Running RelocaTE...")
    # mccutils.run_command(command, log=log)
    print(" ".join(command))
    mccutils.run_command(command)
    combine_gffs(out_dir, out_gff)




def make_annotation_file(te_gff, out_dir):
    annotation_file = out_dir+"/annotation.tsv"
    with open(annotation_file,"w") as out:
        with open(te_gff,"r") as gff:
            for line in gff:
                if "#" not in line:
                    split_line = line.split("\t")
                    outline = "\t".join([split_line[2], split_line[0]+":"+split_line[3]+".."+split_line[4]])
                    out.write(outline+"\n")
    
    return annotation_file


def combine_gffs(out_dir, outgff):
    with open(outgff,"w") as out: 
        for a in os.listdir(out_dir):
            if os.path.isdir(out_dir+"/"+a):
                for b in os.listdir(out_dir+"/"+a):
                    if "results" in b:
                        for c in os.listdir(out_dir+"/"+a+"/"+b):
                            if ".gff" in c:
                                with open(out_dir+"/"+a+"/"+b+"/"+c, "r") as ingff:
                                    for line in ingff:
                                        if "#" not in line:
                                            out.write(line)

if __name__ == "__main__":                
    main()