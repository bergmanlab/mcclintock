import os
import sys
import subprocess
import math
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot
import matplotlib.patches as mpatches
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.coverage.coverage as config


def main():
    mcc_out = snakemake.config["args"]['out']
    mccutils.mkdir(mcc_out+"/results/")
    coverage_out = mcc_out+"/results/coverage/"
    mccutils.mkdir(coverage_out)

    # ensures intermediate files from previous runs are removed
    for f in os.listdir(coverage_out):
        mccutils.remove(coverage_out+"/"+f)

    run_id = snakemake.config['args']['run_id']
    te_seqs = snakemake.input.consensus
    log = snakemake.params.log
    
    # always use consensus fasta for masking the genome
    mccutils.mkdir(coverage_out+"/input")
    mccutils.mkdir(coverage_out+"/te-depth-files")
    masked_reference, masked_gff = repeatmask_genome(snakemake.input.ref, te_seqs, snakemake.threads, run_id, coverage_out, log)

    # uses coverage fasta (if exists) for augmenting and coverage analysis
    if snakemake.config['in']['coverage_fasta'] != "None":
        mccutils.run_command(["cp", snakemake.input.coverage_fa, coverage_out+"/input/coverage_tes_unfixed.fasta"])
        te_seqs = coverage_out+"/input/coverage_tes_unfixed.fasta"
        mccutils.replace_special_chars_fasta(coverage_out+"/input/coverage_tes_unfixed.fasta", te_seqs)

    augmented_reference = augment_genome(masked_reference, te_seqs, run_id, coverage_out)
    index_genome(snakemake.input.ref, log)
    index_genome(augmented_reference, log)
    
    if snakemake.config['in']['fq2'] == "None":
        sam = map_reads(augmented_reference, snakemake.input.fq1, snakemake.threads, snakemake.params.sample, run_id, coverage_out, log)
    else:
        sam = map_reads(augmented_reference, snakemake.input.fq1, snakemake.threads, snakemake.params.sample, run_id, coverage_out, log, fq2=snakemake.input.fq2)

    bam = sam_to_bam(sam, augmented_reference, snakemake.params.sample, snakemake.threads, run_id, coverage_out, log)
    nonte_bed = make_nonte_bed(snakemake.input.ref, masked_gff, run_id, coverage_out, log)
    genome_depth = get_genome_depth(nonte_bed, bam, run_id, coverage_out, log)

    edge_trim = 0
    if config.OMIT_EDGES:
        if config.OMIT_EDGES_READ_LENGTH:
            edge_trim = mccutils.estimate_read_length(snakemake.input.fq1)
        else:
            edge_trim = config.OMIT_EDGES_LENGTH
            
    te_names, all_coverage_files, uniq_coverage_files, avg_norm_te_depths = make_depth_table(te_seqs, bam, genome_depth, run_id, coverage_out, snakemake.output[0], log, trim_edges=edge_trim)
    make_plots(te_names, all_coverage_files, uniq_coverage_files, avg_norm_te_depths, genome_depth, snakemake.params.sample, coverage_out, trim_edges=edge_trim)

    mccutils.remove(sam)
    mccutils.remove(bam)


def repeatmask_genome(reference, lib, threads, run_id, out, log):
    mccutils.log("coverage","Running RepeatMasker",log=log)
    outdir = out+"/input/repeatmasker_"+run_id
    mccutils.mkdir(outdir)
    os.chdir(outdir)
    ref_name = os.path.basename(reference)
    command = ["RepeatMasker","-pa", str(threads), "-lib", lib, "-dir", outdir, "-s", "-gff", "-nolow", "-no_is", reference]
    mccutils.run_command(command, log=log)

    masked_fasta = outdir+"/"+ref_name+".masked"
    repeat_masker_gff = outdir+"/"+ref_name+".out.gff"

    return masked_fasta, repeat_masker_gff


def augment_genome(fasta1, fasta2, run_id, out):
    mccutils.log("coverage","augmenting reference genome")
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

def index_genome(fasta, log):
    mccutils.log("coverage","samtools and bwa indexing reference", log=log)
    mccutils.run_command(["samtools", "faidx", fasta], log=log)
    mccutils.run_command(["bwa", "index", fasta], log=log)


def map_reads(reference, fq1, threads, sample_name, run_id, out, log, fq2=None):
    mccutils.log("coverage","mapping reads to augmented reference genome", log=log)
    command = ["bwa", "mem", "-t", str(threads), "-R", "@RG\\tID:"+sample_name+"\\tSM:"+sample_name, reference, fq1]

    if fq2 is not None:
        command.append(fq2)
    
    sam = out+"/input/"+run_id+"_"+sample_name+".sam"
    mccutils.run_command_stdout(command, sam, log=log)

    return sam

def sam_to_bam(sam, reference, sample_name, threads, run_id, out, log):
    mccutils.log("coverage","converting SAM to BAM, and indexing", log=log)
    threads = str(threads)
    tmp_bam = out+"/input/"+run_id+"_tmp.bam"
    command = ["samtools", "view", "-Sb", "-@", threads, "-t", reference+".fai", sam]
    mccutils.run_command_stdout(command, tmp_bam, log=log)

    sorted_bam = out+"/input/"+run_id+"_"+sample_name+".bam"
    command = ["samtools", "sort", "-@", threads, tmp_bam]
    mccutils.run_command_stdout(command, sorted_bam, log=log)

    mccutils.run_command(["samtools", "index", sorted_bam], log=log)

    mccutils.remove(tmp_bam)

    return sorted_bam

def make_nonte_bed(reference, masked_gff, run_id, out, log):
    mccutils.log("coverage","creating BED file of non-TE regions", log=log)
    masked_bed = out+"/input/"+run_id+"_ref_tes.bed"
    repeatmasker_gff_to_bed(masked_gff, masked_bed)

    sorted_bed = out+"/input/"+run_id+"_ref_tes_sorted.bed"
    mccutils.run_command_stdout(["bedtools", "sort", "-i", masked_bed], sorted_bed, log=log)

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
    mccutils.run_command_stdout(command, non_te_bed, log=log)

    for f in [masked_bed, sorted_bed, chrom_idx]:
        mccutils.remove(f)

    return non_te_bed



def repeatmasker_gff_to_bed(gff, bed):
    mccutils.log("coverage","converting repeatmasker GFF to BED")
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


def get_genome_depth(non_te_bed, bam, run_id, out, log):
    mccutils.log("coverage","determining the coverage depth of the genome", log=log)
    depth_file = out+"/input/"+run_id+"genome.depth"
    command = ["samtools", "depth", "-aa", "-b", non_te_bed, bam, "-d", "0"]
    mccutils.run_command_stdout(command, depth_file, log=log)

    genome_depth = get_avg_depth(depth_file)

    mccutils.remove(depth_file)

    return genome_depth


def get_avg_depth(depth_file, trim_edges=0):
    total = 0
    positions = 0

    depth_lines = []
    with open(depth_file, "r") as depth:
        for line in depth:
            depth_lines.append(line)
    
    for x in range(trim_edges, len(depth_lines)-trim_edges):
            line = depth_lines[x]
            positions += 1
            split_line = line.split("\t")
            total += int(split_line[2])
    
    if positions > 0:
        avg_depth = total/positions
    else:
        mccutils.log("coverage", "ERROR: no positions in: "+depth_file)
        mccutils.log("coverage", "It's likely the length of the family for: "+depth_file+" is shorter than the trimmed edges parameter")
        mccutils.log("coverage", "Setting avg depth to zero for this family, to avoid this, lower the OMIT_EDGES_LENGTH parameter in the coverage config file")
        avg_depth = 0

    return avg_depth

def make_depth_table(te_fasta, bam, genome_depth, run_id, out, depth_csv, log, trim_edges=0):
    mccutils.log("coverage","creating TE depth coverage table", log=log)
    with open(depth_csv, "w") as table:
            table.write("TE-Family,Normalized-Depth,Normalized-Unique-Depth"+"\n")
    
    te_names = []
    uniq_coverage_files = []
    all_coverage_files = []
    avg_norm_depths = []
    avg_uniq_norm_depths = []

    with open(te_fasta,"r") as fa:
        for line in fa:
            if ">" in line:
                te_name = line.replace("\n","")
                te_name = te_name.replace(">","")

                mccutils.mkdir(out+"/te-depth-files")
                highQ = out+"/te-depth-files/"+te_name+".highQ.cov"
                command = ["samtools", "depth", "-aa", "-r", te_name, bam, "-d", "0", "-Q", "1"]
                mccutils.run_command_stdout(command, highQ, log=log)

                allQ = out+"/te-depth-files/"+te_name+".allQ.cov"
                command = ["samtools", "depth", "-aa", "-r", te_name, bam, "-d", "0", "-Q", "0"]
                mccutils.run_command_stdout(command, allQ, log=log)

                # make normalized coverage files
                allQ_chrom, allQ_pos, allQ_cov = read_samtools_depth_file(allQ)
                with open(out+"/te-depth-files/"+te_name+".allQ.normalized.cov","w") as covfile:
                    for i,pos in enumerate(allQ_pos):
                        cov = str(round(allQ_cov[i]/genome_depth,2))
                        line = "\t".join([allQ_chrom,str(pos),cov])
                        covfile.write(line+"\n")
                
                highQ_chrom, highQ_pos, highQ_cov = read_samtools_depth_file(highQ)
                with open(out+"/te-depth-files/"+te_name+".highQ.normalized.cov","w") as covfile:
                    for i,pos in enumerate(highQ_pos):
                        cov = str(round(highQ_cov[i]/genome_depth,2))
                        line = "\t".join([highQ_chrom,str(pos),cov])
                        covfile.write(line+"\n")

                avg_depth = get_avg_depth(allQ, trim_edges=trim_edges)
                avg_norm_depth = avg_depth/genome_depth

                avg_uniq_depth = get_avg_depth(highQ, trim_edges=trim_edges)
                avg_uniq_norm_depth = avg_uniq_depth/genome_depth

                with open(depth_csv, "a") as table:
                    table.write(te_name+","+str(round(avg_norm_depth,2))+","+str(round(avg_uniq_norm_depth,2))+"\n")
    
                te_names.append(te_name)
                uniq_coverage_files.append(highQ)
                all_coverage_files.append(allQ)
                avg_norm_depths.append(avg_norm_depth)
    
    return te_names, all_coverage_files, uniq_coverage_files, avg_norm_depths


def make_plots(te_names, all_coverage_files, uniq_coverage_files, avg_norm_te_depths, genome_depth, sample_name, out, trim_edges=0):
    mccutils.log("coverage","creating TE coverage plots")
    mccutils.mkdir(out+"/plots")
    for x, te_name in enumerate(te_names):
        chrom, all_pos, all_cov = read_samtools_depth_file(all_coverage_files[x])
        chrom2, uniq_pos, uniq_cov = read_samtools_depth_file(uniq_coverage_files[x])

        plot_height = 3
        plot_width = 10
        hline = avg_norm_te_depths[x]
        output = out+"plots/"+te_name+".png"
        plot = plot_coverage(chrom, all_pos, all_cov, uniq_pos, uniq_cov, sample_name, plot_height, plot_width, genome_depth, hline, trim_edges=trim_edges)
        plot.savefig(output, bbox_inches="tight")
        plot.close()
        mccutils.log("coverage","plot created: "+output)



def read_samtools_depth_file(depth_file):
    chrom = ""
    x = []
    y = []
    
    with open(depth_file,"r") as depth:
        for line in depth:
            split_line = line.split("\t")
            if chrom == "":
                chrom = split_line[0]
            elif chrom != split_line[0]:
                sys.exit("Error: Samtools depth file has multiple reference contigs.... exiting...\n")
            
            x.append(int(split_line[1]))
            y.append(int(split_line[2]))

    return chrom, x, y


def plot_coverage(chrom, all_pos, all_cov, uniq_pos, uniq_cov, title, height, width, normalize_cov, add_hline, trim_edges=0 ):

    if normalize_cov:
        all_cov = np.array(all_cov)
        uniq_cov = np.array(uniq_cov)

        all_cov = all_cov/normalize_cov
        uniq_cov = uniq_cov/normalize_cov

    fig = matplotlib.pyplot.figure(figsize=(width, height), dpi=300)

    plot = matplotlib.pyplot

    ax = plot.subplot(1,1,1)

    ax.fill_between(all_pos, all_cov, np.zeros(len(all_cov)), color='grey', step="pre", alpha=0.15)

    ax.fill_between(uniq_pos, uniq_cov, np.zeros(len(uniq_cov)), color='darkgrey', step="pre", alpha=0.4)


    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(axis='y', colors='grey')
    ax.tick_params(axis='x', labelsize=5, length=0, pad=5)
    ax.tick_params(axis='y', labelsize=5, length=0, pad=5)

    x_min = min(all_pos)
    x_max = max(all_pos)

    y_min = min(all_cov)
    y_max = max(all_cov)

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(0, y_max)
    plot.xticks([x_min, int(np.median(all_pos)), x_max])

    y_max_tick = y_max
    y_max_tick = math.ceil(y_max_tick*100)/100


    if add_hline is None:
        y_mid_tick = y_max/2
        y_mid_tick = math.ceil(y_mid_tick*100)/100
    else:
        y_mid_tick = add_hline

    plot.yticks([0, y_mid_tick, y_max_tick])
    
    ax.set_xlabel('Position on '+ chrom, fontsize=8)

    # prevents scientific notation if positions are large
    ax.ticklabel_format(useOffset=False, style='plain')
    

    if title:
        ax.set_title(" "*15 + title, fontsize=10, loc='left',y=1.07)


    ## legend

    labels = []
    elements = []

    elements += [mpatches.Patch(color='darkgrey',alpha=.4)]
    labels.append("Coverage (MAPQ > 0)")
    elements += [mpatches.Patch(color='grey',alpha=.15)]
    labels.append("Coverage (MAPQ = 0)")

    if add_hline is not None:
        ax.set_ylabel('Normalized Coverage', fontsize=8)
        ax.axhline(xmin=(x_min+trim_edges)/x_max, xmax=(x_max - trim_edges)/x_max, y=add_hline, color='black', alpha=0.5)
        elements += [matplotlib.pyplot.Line2D([0,0],[0,1],color='black', alpha=0.5)]
        labels.append("Avg. Norm. Cov.")
    else:
        ax.set_ylabel('Coverage', fontsize=8)

    plot.legend( elements , labels, loc=9, fontsize = 6, frameon=False, ncol=3, bbox_to_anchor=(0.5, 1.1))


    return plot

if __name__ == "__main__":                
    main()