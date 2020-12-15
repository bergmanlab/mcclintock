import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    mccutils.log("teflon","setting up for TEFLoN")

    te_gff = snakemake.input.te_gff
    taxonomy = snakemake.input.taxonomy
    consensus = snakemake.input.consensus
    reference_genome = snakemake.input.reference_genome
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2

    threads = snakemake.threads
    out_dir = snakemake.params.out_dir
    script_dir = snakemake.params.script_dir
    log = snakemake.params.log

    ref_bed = snakemake.output.ref_bed
    teflon_taxonomy = snakemake.output.teflon_taxonomy

    make_reference_bed(te_gff, ref_bed)

    make_taxonomy_file(taxonomy, teflon_taxonomy)

    prep_annotations(script_dir, out_dir, ref_bed, teflon_taxonomy, consensus, reference_genome, log=log)

    map_reads(out_dir, fq1, fq2, threads=threads, log=log)

    mccutils.log("teflon","setup for TEFLoN complete")



def make_reference_bed(reference_gff, reference_bed):
    bed_lines = []
    with open(reference_gff, "r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                contig = split_line[0]
                start = str(int(split_line[3])-1)
                end = split_line[4]
                strand = split_line[6]
                feats = split_line[8].replace("\n","")
                te_id = ""
                split_feats = feats.split(";")
                for feat in split_feats:
                    if "ID=" == feat[0:3]:
                        te_id = feat.replace("ID=","")

                bed_lines.append("\t".join([contig, start, end, te_id, ".", strand, "."]) + "\n")

    with open(reference_bed, "w") as bed:
        for line in bed_lines:
            bed.write(line)


def make_taxonomy_file(taxonomy, teflon_taxonomy):
    taxonomy_lines = ["id\tfamily\n"]

    with open(taxonomy, "r") as tsv:
        for line in tsv:
            taxonomy_lines.append(line)
    
    with open(teflon_taxonomy, "w") as out:
        for line in taxonomy_lines:
            out.write(line)


def prep_annotations(script_dir, out_dir, ref_bed, taxonomy, consensus, reference, log=None):
    command = [
        "python", script_dir+"teflon_prep_annotation.py",
        "-wd", out_dir,
        "-a", ref_bed,
        "-t", taxonomy,
        "-f", consensus,
        "-g", reference,
        "-p", "teflon"
    ]

    mccutils.run_command(command, log=log)


def map_reads(out_dir, fq1, fq2, threads=1, log=None):
    reference_genome = out_dir+"/teflon.prep_MP/teflon.mappingRef.fa"
    command = ["bwa", "index", reference_genome]
    mccutils.run_command(command, log=log)

    out_sam = out_dir+"teflon.sam"

    command = [
        "bwa", "mem",
        "-t", str(threads),
        "-Y", reference_genome,
        fq1, 
        fq2
    ]

    mccutils.run_command_stdout(command, out_sam, log=log)

    out_bam = out_dir+"teflon.bam"
    command = ["samtools", "view", "-Sb", out_sam]
    mccutils.run_command_stdout(command, out_bam, log=log)

    sorted_bam = out_dir+"teflon.sorted.bam"
    command = ["samtools", "sort", "-@", str(threads), "-o", sorted_bam, out_bam]
    mccutils.run_command(command, log=log)

    command = ["samtools", "index", sorted_bam ]
    mccutils.run_command(command, log=log)

    mccutils.remove(out_sam)
    mccutils.remove(out_bam)

    return sorted_bam



if __name__ == "__main__":                
    main()