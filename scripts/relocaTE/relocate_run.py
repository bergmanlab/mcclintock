import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.relocate.relocate_run as config


def main():
    

    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    raw_fq2 = snakemake.params.raw_fq2
    is_paired = True
    if raw_fq2 == "None":
        is_paired = False

    script_dir = snakemake.params.script_dir
    out_dir = snakemake.params.out_dir
    out_gff = snakemake.output[0]

    mccutils.log("relocate","running RelocaTE", log=log)

    input_dir = snakemake.params.out_dir+"/input/"
    mccutils.remove(input_dir)
    mccutils.mkdir(input_dir)
    fq_dir = input_dir+"fastq/"
    mccutils.mkdir(fq_dir)

    consensus_fasta = input_dir+"consensus.fasta"
    te_gff = input_dir+"te.gff"
    reference_fasta = input_dir+"reference.fasta"

    os.symlink(snakemake.input.consensus_fasta, consensus_fasta)
    os.symlink(snakemake.input.te_gff, te_gff)
    os.symlink(snakemake.input.reference_fasta, reference_fasta)
    if is_paired:
        os.symlink(snakemake.input.fq1, fq_dir+sample_name+"_1.fq")
        os.symlink(snakemake.input.fq2, fq_dir+sample_name+"_2.fq")
    else:
        os.symlink(snakemake.input.fq1, fq_dir+sample_name+".unPaired.fq")




    annotation = make_annotation_file(te_gff, out_dir)
    os.chdir(out_dir)

    command = ["perl", script_dir+"/relocaTE.pl", 
                    "-t", consensus_fasta, 
                    "-d", fq_dir, 
                    "-g", reference_fasta, 
                    "-o", ".", 
                    "-r", annotation,
                    "-l", str(config.RELOCATE['l']),
                    "-m", str(config.RELOCATE['m']),
                    "-bm", str(config.RELOCATE['bm']),
                    "-bt", str(config.RELOCATE['bt']),
                    "-f", str(config.RELOCATE['f'])]


    if is_paired:
        command += ["-1", "_1", "-2", "_2"]
    else:
        command += ["-u", "unPaired"]
    
    
    mccutils.run_command(command, log=log)
    combine_gffs(out_dir, out_gff)
    mccutils.log("relocate","RelocaTE run complete")




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