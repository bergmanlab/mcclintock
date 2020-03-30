import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils


def main():
    mcc_out = snakemake.config["args"]['out']
    coverage_out = snakemake.config["args"]['out']+"/results/coverage/"
    run_id = snakemake.config['args']['run_id']
    te_seqs = snakemake.input.consensus
    if snakemake.config['in']['coverage_fasta'] != "None":
        te_seqs = snakemake.input.coverage_fa
    
    mccutils.mkdir(coverage_out+"/input")
    masked_reference, masked_gff = repeatmask_genome(snakemake.input.ref, te_seqs, snakemake.threads, run_id, coverage_out)
    augmented_reference = augment_genome(masked_reference, te_seqs, run_id, coverage_out)
    index_genome(snakemake.input.ref)
    index_genome(augmented_reference)
    
    if snakemake.config['in']['fq2'] != "None":
        sam = map_reads(augmented_reference, snakemake.input.fq1, snakemake.threads, snakemake.params.sample, run_id, coverage_out)
    else:
        sam = map_reads(augmented_reference, snakemake.input.fq1, snakemake.threads, snakemake.params.sample, run_id, coverage_out, fq2=snakemake.input.fq2)

    bam = sam_to_bam(sam, augmented_reference, snakemake.params.sample, snakemake.threads, run_id, coverage_out)

    nonte_bed = make_nonte_bed(snakemake.input.ref, augmented_reference, masked_gff, run_id, coverage_out)

    genome_depth = get_genome_depth(nonte_bed, bam, run_id, coverage_out)

    make_depth_table(te_seqs, bam, genome_depth, run_id, coverage_out)

    mccutils.run_command(["touch", snakemake.config['args']['out']+"/coverage/coverage.log"])

def repeatmask_genome(reference, lib, threads, run_id, out):
    outdir = out+"/input/repeatmasker_"+run_id
    mccutils.mkdir(outdir)
    os.chdir(outdir)

    command = ["RepeatMasker","-pa", str(threads), "-lib", lib, "-dir", outdir, "-s", "-gff", "-nolow", "-no_is", reference]
    mccutils.run_command(command)

    ref_name = os.path.basename(reference)
    masked_fasta = outdir+"/"+ref_name+".masked"
    repeat_masker_gff = outdir+"/"+ref_name+".out.gff"

    return masked_fasta, repeat_masker_gff


def augment_genome(fasta1, fasta2, run_id, out):
    augmented_genome = out+"/input/"+run_id+"_augmented_reference.fasta"
    lines = []
    with open(fasta1,"r") as fa1:
        for line in fa1:
            lines.append(line)
    
    with open(fasta2, "r") as fa2:
        for line in fa2:
            lines.append(line)
    
    with open(augmented_genome, "w") as outfa:
        for line in lines:
            outfa.write(line)
    
    return augmented_genome

def index_genome(fasta):
    mccutils.run_command(["samtools", "faidx", fasta])
    mccutils.run_command(["bwa", "index", fasta])


def map_reads(reference, fq1, threads, sample_name, run_id, out, fq2=None):
    command = ["bwa", "mem", "-t", str(threads), "-R", "@RG\\tID:"+sample_name+"\\tSM:"+sample_name, reference, fq1]

    if fq2 is not None:
        command.append(fq2)
    
    sam = out+"/input/"+run_id+"_"+sample_name+".sam"
    mccutils.run_command_stdout(command, sam)

    return sam

def sam_to_bam(sam, reference, sample_name, threads, run_id, out):
    threads = str(threads)
    tmp_bam = out+"/input/"+run_id+"_tmp.bam"
    command = ["samtools", "view", "-Sb", "-@", threads, "-t", reference+".fai", sam]
    mccutils.run_command_stdout(command, tmp_bam)

    sorted_bam = out+"/input/"+run_id+"_"+sample_name+".bam"
    command = ["samtools", "sort", "-@", threads, tmp_bam]
    mccutils.run_command_stdout(command, sorted_bam)

    mccutils.run_command(["samtools", "index", sorted_bam])

    mccutils.run_command(["rm", tmp_bam])

    return sorted_bam

def make_nonte_bed(reference, augmented_reference, masked_gff, run_id, out):
    masked_bed = out+"/input/"+run_id+"_ref_tes.bed"
    repeatmasker_gff_to_bed(masked_gff, masked_bed)

    sorted_bed = out+"/input/"+run_id+"_ref_tes_sorted.bed"
    mccutils.run_command_stdout(["bedtools", "sort", "-i", masked_bed], sorted_bed)

    chromosome_names = []
    with open(reference, "r") as fa:
        for line in fa:
            if ">" in line:
                chromosome_names.append(line.replace(">","").replace("\n",""))

    chrom_idx = out+"/input/"+run_id+"_ref.genome"
    with open(reference+".fai", "r") as faidx:
        with open(chrom_idx,"w") as genome:
            for line in faidx:
                split_line = line.split("\t")
                out_line = "\t".join([split_line[0],split_line[1]])
                genome.write(out_line+"\n")
    
    non_te_bed = out+"/input/"+run_id+"_ref_nonte.bed"
    command = ["bedtools", "complement", "-i", sorted_bed, "-g", chrom_idx]
    mccutils.run_command_stdout(command, non_te_bed)


    mccutils.run_command(["rm", masked_bed, sorted_bed, chrom_idx])

    return non_te_bed



def repeatmasker_gff_to_bed(gff, bed):
    with open(gff,"r") as ingff:
        with open(bed, "w") as outbed:
            for line in ingff:
                if "#" not in line:
                    split_line = line.split("\t")
                    feats = split_line[8].split(" ")
                    te = "MISSING"
                    for feat in feats:
                        if "Motif" in feat:
                            te = feat.split(":")[1]
                            te = te.replace('"','')
                    out_line = "\t".join([split_line[0], str(int(split_line[3])-1), split_line[4], te, ".", split_line[6]])
                    outbed.write(out_line+"\n")


def get_genome_depth(non_te_bed, bam, run_id, out):
    depth_file = out+"/input/"+run_id+"genome.depth"
    command = ["samtools", "depth", "-aa", "-b", non_te_bed, bam, "-d", "0"]
    mccutils.run_command_stdout(command, depth_file)

    genome_depth = get_avg_depth(depth_file)

    mccutils.run_command(["rm", depth_file])

    return genome_depth


def get_avg_depth(depth_file):
    total = 0
    positions = 0

    with open(depth_file, "r") as depth:
        for line in depth:
            positions += 1
            split_line = line.split("\t")

            total += int(split_line[2])
    
    avg_depth = total/positions

    return avg_depth

def make_depth_table(te_fasta, bam, genome_depth, run_id, out):
    depth_csv = out+"/output/te_depth.csv"
    with open(depth_csv, "w") as table:
            table.write("TE-Family"+"\t"+"Normalized-Depth"+"\n")

    with open(te_fasta,"r") as fa:
        for line in fa:
            if ">" in line:
                te_name = line.replace("\n","")
                te_name = te_name.replace(">","")

                mccutils.mkdir(out+"/output/te-depth")
                highQ = out+"/output/te-depth/"+te_name+".highQ.cov"
                command = ["samtools", "depth", "-aa", "-r", te_name, bam, "-d", "0", "-Q", "1"]
                mccutils.run_command_stdout(command, highQ)

                allQ = out+"/output/te-depth/"+te_name+".allQ.cov"
                command = ["samtools", "depth", "-aa", "-r", te_name, bam, "-d", "0", "-Q", "0"]
                mccutils.run_command_stdout(command, allQ)

                avg_depth = get_avg_depth(allQ)
                avg_norm_depth = avg_depth/genome_depth

                with open(depth_csv, "a") as table:
                    table.write(te_name+"\t"+str(avg_norm_depth)+"\n")
    


if __name__ == "__main__":                
    main()