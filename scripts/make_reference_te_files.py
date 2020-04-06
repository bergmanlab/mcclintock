#!/usr/bin/env python3

import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import modules.mccutils as mccutils
import modules.fix_fasta as fix_fasta
from Bio import SeqIO


def main():
    reference_fa = snakemake.input[0]
    consensus_TEs = snakemake.input[1]
    locations_gff = snakemake.input[2]
    taxonomy_tsv = snakemake.input[3]
    ref_tes = snakemake.config['in']['locations']
    processors = snakemake.config['args']['proc']
    mcc_out = snakemake.config['args']['out']
    run_id = snakemake.config['args']['run_id']

    print("<PROCESSING> creating reference TE files...")
    # if no reference TEs were provided, they must be found with repeatmasker
    if ref_tes == "None":
        print("<PROCESSING> no reference TEs provided... finding reference TEs with RepeatMasker...")
        masked_reference, formatted_ref_tes = repeat_mask(reference_fa, consensus_TEs, processors, run_id, mcc_out)
        formatted_consensus_TEs, te_families = format_consensus_tes(consensus_TEs, run_id, mcc_out)
        formatted_taxonomy_tsv, formatted_ref_tes = make_te_taxonomy_map(formatted_ref_tes, te_families, run_id, mcc_out)

    # use provided reference TEs to mask genome
    else:
        formatted_ref_tes = format_ref_te_gff(ref_tes, run_id, mcc_out)
        masked_reference = mask_reference(reference_fa, formatted_ref_tes, run_id, mcc_out)
        formatted_taxonomy_tsv = taxonomy_tsv
        formatted_consensus_TEs = consensus_TEs


    ref_te_fasta = get_ref_te_fasta(reference_fa, formatted_ref_tes, run_id, mcc_out)

    popoolationTE_tsv = make_popoolationTE_taxonomy(formatted_taxonomy_tsv, formatted_consensus_TEs, run_id, mcc_out)

    mccutils.run_command(["cp", formatted_ref_tes, snakemake.output[7]])

    augmented_reference = augment_reference(reference_fa, ref_te_fasta, formatted_consensus_TEs, run_id, mcc_out, 
                                            add_consensus=eval(snakemake.config['args']['include_consensus']), 
                                            add_reference=eval(snakemake.config['args']['include_reference']))
    
    formatted_ref_tes = augment_ref_te_gff(formatted_ref_tes, ref_te_fasta, formatted_consensus_TEs, run_id, mcc_out,
                                            add_consensus=eval(snakemake.config['args']['include_consensus']),
                                            add_reference=eval(snakemake.config['args']['include_reference']))
    
    formatted_taxonomy_tsv = augment_taxonomy_map(formatted_taxonomy_tsv, ref_te_fasta, formatted_consensus_TEs, run_id, mcc_out,
                                            add_consensus=eval(snakemake.config['args']['include_consensus']),
                                            add_reference=eval(snakemake.config['args']['include_reference']))


    mccutils.run_command(["mv", masked_reference, snakemake.output[0]])
    mccutils.run_command(["mv", formatted_ref_tes, snakemake.output[1]])
    mccutils.run_command(["mv", formatted_taxonomy_tsv, snakemake.output[2]])
    mccutils.run_command(["cp", formatted_consensus_TEs, snakemake.output[3]])
    mccutils.run_command(["mv", ref_te_fasta, snakemake.output[4]])
    mccutils.run_command(["mv", augmented_reference, snakemake.output[5]])
    mccutils.run_command(["mv", popoolationTE_tsv, snakemake.output[6]])

    print("<PROCESSING> reference TE files created")



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
                if "#" not in line:
                    line = line.replace("McClintock-int","McClintock")
                    line = line.replace("POGON1", "pogo")
                    split_line = line.split("\t")
                    feats = split_line[8]
                    te = feats.split(" ")[1]
                    te = te.replace('"','').split(":")[1]
                    feats = ";".join(["ID="+te, "Name="+te, "Alias="+te])
                    split_line[2] = te
                    split_line[8] = feats
                    line = "\t".join(split_line)

                    outgff.write(line+'\n')


    masked_fasta = outdir+"/"+ref_name+".masked"
    fasta_lines = fix_fasta.fix_fasta_lines(masked_fasta, 80)

    masked_ref = out+"/tmp/"+run_id+"tmpmaskedreference.fasta"
    with open(masked_ref,"w") as fa:
        for line in fasta_lines:
            fa.write(line+"\n")
    
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


def make_te_taxonomy_map(ref_tes, consensus_families, run_id, out):
    te_family_counts = {}
    family_tsv = out+"/tmp/"+run_id+"tmpfam.tsv"
    gff_lines = []
    with open(ref_tes, "r") as ingff:
        with open(family_tsv,"w") as outtsv:
            for line in ingff:
                if "#" not in line:
                    split_line = line.split("\t")
                    family = split_line[2]
                    if family in te_family_counts.keys():
                        te_family_counts[family] += 1
                    else:
                        te_family_counts[family] = 1
                    
                    te_id = family+"_"+str(te_family_counts[family])
                    split_line[2] = te_id
                    split_line[8] = ";".join(["ID="+te_id, "Name="+te_id, "Alias="+te_id])
                    line = "\t".join(split_line)
                    gff_lines.append(line+"\n")
                    if family in consensus_families:
                        outtsv.write(te_id+"\t"+family+"\n")
    
    with open(ref_tes, "w") as outgff:
        for line in gff_lines:
            outgff.write(line)

    return family_tsv, ref_tes




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

def get_ref_te_fasta(reference, ref_tes_gff, run_id, out):
    log = out+"/logs/bedtools.getfasta.log"
    ref_te_fasta = out+"/tmp/"+run_id+"tmpreferencetes.fasta"
    command = ["bedtools", "getfasta", "-name", "-fi", reference, "-bed", ref_tes_gff, "-fo", ref_te_fasta]
    mccutils.run_command(command, log=log)

    # fasta_lines = fixfasta.fix_fasta_lines(ref_te_fasta, 80)
    # with open(ref_te_fasta, "w") as outfa:
    #     for line in fasta_lines:
    #         outfa.write(line+"\n")

    return ref_te_fasta

def make_popoolationTE_taxonomy(taxonomy, consensus, run_id, out):
    te_families = []
    popoolationTE_taxonomy = out+"/tmp/"+run_id+"tmppopoolationtetaxonomy.tsv"

    with open(taxonomy, "r") as intsv:
        for line in intsv:
            split_line = line.split("\t")
            te_families.append([split_line[0], split_line[1].replace("\n","")])
    
    fasta_records = SeqIO.parse(consensus,"fasta")

    for record in fasta_records:
        element = str(record.id)
        te_families.append([element, element])
    
    with open(popoolationTE_taxonomy, "w") as outtsv:
        header = "\t".join(["insert","id","family","superfamily","suborder","order","class","problem\n"])
        outtsv.write(header)

        for te_fam in te_families:
            element = te_fam[0]
            family = te_fam[1]
            line = "\t".join([element, family, family, family, "na", "na", "na", "0\n"])
            outtsv.write(line)
    
    return popoolationTE_taxonomy



def augment_reference(reference, ref_tes, consensus_tes, run_id, out, add_consensus=False, add_reference=False):
    out_lines = fix_fasta.fix_fasta_lines(reference, 80)
    if add_consensus:
        consensus_lines = fix_fasta.fix_fasta_lines(consensus_tes, 80)
        out_lines += consensus_lines
    if add_reference:
        reference_te_lines = fix_fasta.fix_fasta_lines(ref_tes, 80)
        out_lines += reference_te_lines
    
    augmented_reference = out+"/tmp/"+run_id+"tmpaugmentedreference.fasta"
    with open(augmented_reference,"w") as outfasta:
        for line in out_lines:
            outfasta.write(line+"\n")
    
    return augmented_reference


def augment_ref_te_gff(ref_tes_gff, ref_te_fasta, consensus_te_fasta, run_id, out, add_consensus=False, add_reference=False):
    gff_lines = []

    with open(ref_tes_gff,"r") as ingff:
        for line in ingff:
            gff_lines.append(line)
    
    if add_consensus:
        fasta_records = SeqIO.parse(consensus_te_fasta,"fasta")
        for record in fasta_records:
            length = len(str(record.seq))
            te = str(record.id)
            features = ";".join(["ID=instance"+te,"Name=instance"+te,"Alias=instance"+te])
            line = "\t".join([te, "reannotate", "transposable_element", "1", str(length), ".", "+", ".", features])
            line = line+"\n"
            gff_lines.append(line)
    
    if add_reference:
        fasta_records = SeqIO.parse(ref_te_fasta,"fasta")
        for record in fasta_records:
            length = len(str(record.seq))
            te = str(record.id)
            features = ";".join(["ID=instance"+te,"Name=instance"+te,"Alias=instance"+te])
            line = "\t".join([te, "reannotate", "transposable_element", "1", str(length), ".", "+", ".", features])
            line = line+"\n"
            gff_lines.append(line)
    
    with open(ref_tes_gff, "w") as outgff:
        for line in gff_lines:
            outgff.write(line)
    
    return ref_tes_gff


def augment_taxonomy_map(family_map, ref_te_fasta, consensus_te_fasta, run_id, out, add_consensus=False, add_reference=False):
    families = {}
    map_lines = []
    with open(family_map, "r") as tsv:
        for line in tsv:
            map_lines.append(line)
            split_line = line.split("\t")
            families[split_line[0]] = split_line[1].replace("\n", "")
    

    if add_consensus:
        fasta_records = SeqIO.parse(consensus_te_fasta,"fasta")
        for record in fasta_records:
            te = str(record.id)
            line = te+"\t"+te+"\n"
            map_lines.append(line)
    
    # if add_reference:
    #     fasta_records = SeqIO.parse(ref_te_fasta,"fasta")
    #     for record in fasta_records:
    #         te = str(record.id)
    #         family = families[te]
    #         line = te+"\t"+family+"\n"
    #         map_lines.append(line)

    with open(family_map, "w") as outtsv:
        for line in map_lines:
            outtsv.write(line)

    return family_map    


if __name__ == "__main__":                
    main()