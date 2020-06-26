import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.retroseq.retroseq_run as config
import scripts.fix_fasta as fix_fasta
from Bio import SeqIO


def main():
    consensus_fasta = snakemake.input.consensus_fasta
    bam = snakemake.input.bam
    ref_fasta = snakemake.input.ref_fasta
    ref_te_bed = snakemake.input.ref_te_bed
    taxonomy = snakemake.input.taxonomy
    log = snakemake.params.log

    with open(log,"a") as l:
        l.write("consensus fasta: "+consensus_fasta+"\n")
        l.write("BAM: "+bam+"\n")
        l.write("reference fasta: "+ref_fasta+"\n")
        l.write("taxonomy TSV: "+ taxonomy+"\n")
        

    script_dir = snakemake.params.script_dir
    out_dir = snakemake.params.out_dir
    ref_name = snakemake.params.ref_name
    sample_name = snakemake.params.sample_name

    mccutils.log("retroseq","running RetroSeq", log=log)    

    elements = split_consensus_fasta(consensus_fasta, ref_name, out_dir)

    bed_location_file = make_consensus_beds(elements, ref_name, ref_te_bed, taxonomy, out_dir)

    run_retroseq(bam, bed_location_file, ref_fasta, script_dir, sample_name, out_dir, config.PARAMETERS, log=log)
    mccutils.log("retroseq","RetroSeq complete")



def split_consensus_fasta(fasta, ref_name, out):
    elements = []
    out_dir = out+"/split_fasta/"
    mccutils.mkdir(out_dir)
    fasta_records = SeqIO.parse(fasta,"fasta")
    for record in fasta_records:
        fasta_name = str(record.id)
        elements.append(fasta_name)
        special_chars = [";","&","(",")","|","*","?","[","]","~","{","}","<","!","^",'"',"'","\\","$","/"]
        for char in special_chars:
            fasta_name = fasta_name.replace(char,"_")

        tmp_fasta = out_dir+ref_name+"_"+fasta_name+".fasta.tmp"
        with open(tmp_fasta,"w") as outfa:
            outfa.write(">"+str(record.id)+"\n")
            outfa.write(str(record.seq)+"\n")
        
        fasta_lines = fix_fasta.fix_fasta_lines(tmp_fasta, 80)
        out_fasta = out_dir+ref_name+"_"+fasta_name+".fasta"
        with open(out_fasta,"w") as outfa:
            for line in fasta_lines:
                outfa.write(line+"\n")
        
        mccutils.remove(tmp_fasta)
    
    return elements


def make_consensus_beds(elements, ref_name, te_bed, taxon, out):
    out_dir = out+"/split_bed/"
    mccutils.mkdir(out_dir)
    taxon_map = {}
    location_file = out_dir+ref_name+".locationlist"
    with open(taxon,"r") as t:
        for line in t:
            split_line = line.split("\t")
            element = split_line[0]
            element_fam = split_line[1].replace("\n","")
            if element_fam in taxon_map.keys():
                taxon_map[element_fam].append(element)
            else:
                taxon_map[element_fam] = [element]
    
    with open(location_file,"w") as locations:
        for fam in taxon_map.keys():
            if fam in elements:
                bed_name = fam
                special_chars = [";","&","(",")","|","*","?","[","]","~","{","}","<","!","^",'"',"'","\\","$","/"]
                for char in special_chars:
                    bed_name = bed_name.replace(char,"_")
                
                bed_name = out_dir+ref_name+"_"+bed_name+".bed"
                locations.write(fam+"\t"+bed_name+"\n")
                with open(bed_name, "w") as outbed:
                    with open(te_bed,"r") as inbed:
                        for line in inbed:
                            split_line = line.split("\t")
                            element_name = split_line[3]
                            if element_name in taxon_map[fam]:
                                outbed.write(line)
            
    
    return location_file
    

def run_retroseq(bam, bed_locations, ref_fasta, script_dir, sample_name, out_dir, params, log=None):
    discovery_out = out_dir+"/"+sample_name+".discovery"
    command = ["perl", script_dir+"/retroseq.pl", "-discover",  "-bam", bam, "-refTEs", bed_locations, "-output", discovery_out, "-depth", str(params["depth"]), "-reads", str(params['reads']), "-q", str(params['q'])]
    mccutils.run_command(command, log=log)

    call_out  = out_dir+"/"+sample_name+".call"
    command = ["perl", script_dir+"/retroseq.pl", "-call", "-bam", bam, "-input", discovery_out, "-filter", bed_locations, "-ref", ref_fasta, "-output", call_out, "-orientate", "yes", "-depth", str(params["depth"]), "-reads", str(params['reads']), "-q", str(params['q'])]
    mccutils.run_command(command, log=log)
    

            

if __name__ == "__main__":                
    main()