#!/usr/bin/env python3

import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils
import scripts.fix_fasta_lines as fixfasta


def main():
    reference_fa = snakemake.input[0]
    consensus_TEs = snakemake.input[1]
    locations_gff = snakemake.input[2]
    family_tsv = snakemake.input[3]
    ref_tes = snakemake.config['in']['locations']
    processors = snakemake.config['args']['proc']
    mcc_out = snakemake.config['args']['out']
    run_id = snakemake.config['args']['run_id']

    # if no reference TEs were provided, they must be found with repeatmasker
    if ref_tes == "None":
        masked_reference, formatted_ref_tes = repeat_mask(reference_fa, consensus_TEs, processors, run_id, mcc_out)
        formatted_consensus_TEs, te_families = format_consensus_tes(consensus_TEs, run_id, mcc_out)
        formatted_family_tsv = make_te_family_map(formatted_ref_tes, te_families, run_id, mcc_out)

    # use provided reference TEs to mask genome
    else:
        formatted_ref_tes = format_ref_te_gff(ref_tes, run_id, mcc_out)
        masked_reference = mask_reference(reference_fa, formatted_ref_tes, run_id, mcc_out)
        formatted_family_tsv = family_tsv
        formatted_consensus_TEs = consensus_TEs

    mccutils.run_command(["mv", masked_reference, snakemake.output[0]])
    mccutils.run_command(["mv", formatted_ref_tes, snakemake.output[1]])
    mccutils.run_command(["mv", formatted_family_tsv, snakemake.output[2]])
    mccutils.run_command(["mv", formatted_consensus_TEs, snakemake.output[3]])

# creates reference TE gff using repeatmasker and consensus TE fasta
def repeat_mask(reference, te_fasta, procs, run_id, out):
    outdir = out+"/tmp/repeatmasker_"+run_id
    mccutils.mkdir(outdir)
    os.chdir(outdir)
    command = ["RepeatMasker","-pa", procs, "-lib", te_fasta, "-dir", outdir, "-s", "-gff", "-nolow", "-no_is", reference]
    log = out+"/logs/repeatmasker.log"
    mccutils.run_command(command, log=log)
    os.chdir(out)

    # RepeatMasker appears to override the custom database names during the ProcessRepeats
    # this step changes them back, more rules may be needed for other reference genomes
    ref_name = os.path.basename(reference)
    repeatmasker_gff = outdir+"/"+ref_name+".out.gff"
    formatted_ref_tes = out+"/tmp/"+run_id+"tmpreferenceTEs.gff"
    with open(repeatmasker_gff,"r") as rmgff:
        with open(formatted_ref_tes,"w") as outgff:
            for line in rmgff:
                line = line.replace("McClintock-int","McClintock")
                line = line.replace("POGON1", "pogo")
                outgff.write(line)


    masked_fasta = outdir+"/"+ref_name+".masked"
    fasta_lines = fixfasta.fix_fasta_lines(masked_fasta, 80)

    masked_ref = out+"/tmp/"+run_id+"tmpmaskedreference.fasta"
    with open(masked_ref,"w") as fa:
        for line in fasta_lines:
            fa.write(line+"\n")
    
    # mccutils.run_command(["rm","-r", outdir])
    return masked_ref, formatted_ref_tes

def format_consensus_tes(consensus_TEs, run_id, out):
    families = []
    formatted_consensus = out+"/tmp/"+run_id+"tmpformattedconsensus.fasta"
    with open(consensus_TEs, "r") as infa:
        with open(formatted_consensus,"w") as outfa:
            for line in infa:
                if ">" in line:
                    line = line.replace("\n","")
                    family = line.replace(">","")
                    family = family.split("#")[0]
                    if family not in families:
                        families.append(family)
                    
                    outfa.write(">"+family+"\n")
                else:
                    outfa.write(line)
    
    return formatted_consensus, families


def make_te_family_map(ref_tes, consensus_families, run_id, out):
    te_family_counts = {}
    family_tsv = out+"/tmp/"+run_id+"tmpfam.tsv"
    with open(ref_tes, "r") as ingff:
        with open(family_tsv,"w") as outtsv:
            for line in ingff:
                if "#" not in line:
                    split_line = line.split("\t")
                    split_feats = split_line[8].split(" ")
                    target = split_feats[1].replace('"','')
                    family = target.split(":")[1]
                    if family in te_family_counts.keys():
                        te_family_counts[family] += 1
                    else:
                        te_family_counts[family] = 1
                    
                    te_id = family+"_"+str(te_family_counts[family])
                    if family in consensus_families:
                        outtsv.write(te_id+"\t"+family+"\n")
    
    return family_tsv




# standardizes format of reference te gff
def format_ref_te_gff(ref_tes, run_id, out):
    format_ref_tes = out+"/tmp/"+run_id+"tmpreferenceTEs.gff"
    with open(ref_tes,"r") as ingff:
        with open(format_ref_tes,"w") as outgff:
            for line in ingff:
                if "#" not in line:
                    split_line = line.split("\t")
                    features = split_line[8].replace("\n","")
                    split_feats = features.split(";")
                    te_id = "MISSING"
                    for feat in split_feats:
                        if "ID=" in feat:
                            te_id = feat.split("=")[1]
                    
                    split_line[2] = te_id
                    features = ";".join(["ID="+te_id, "Name="+te_id, "Alias="+te_id])
                    line = "\t".join(split_line[0:8])
                    line = line+"\t"+features+"\n"
                    outgff.write(line)
    
    return format_ref_tes
    
    
    
# masks reference genome using reference TEs provided by user
def mask_reference(reference, ref_tes_gff, run_id, out):
    masked_reference = out+"/tmp/"+run_id+"tmpmaskedreference.fasta"
    log = out+"/logs/bedtools.maskfasta.log"
    command = ["bedtools", "maskfasta", "-fi", reference, "-fo", masked_reference, "-bed", ref_tes_gff]
    mccutils.run_command(command, log=log)

    return masked_reference

    

if __name__ == "__main__":                
    main()