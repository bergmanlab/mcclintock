import argparse
import os
import sys
import traceback
import json
import random
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict

def main():
    args = parse_args()
    consensus_seqs = get_seqs(args.consensus)
    reference_seqs = get_seqs(args.reference)
    for x in range(0,args.reps):
        modified_reference = add_synthetic_insertion(reference_seqs, consensus_seqs, args.config, x, args.out, run_id=args.runid, seed=args.seed)

def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock Simulation', description="Script to run synthetic insertion simulations to evaluate mcclintock component methods")

    ## required ##
    parser.add_argument("--mode", type=str, help="which mode this script should be run", required=True)
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required='run' in sys.argv)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required='run' in sys.argv)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required='run' in sys.argv)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required='run' in sys.argv)
    parser.add_argument("-j","--config", type=str, help="A json config file containing information on TE family TSD size and target sites", required='run' in sys.argv)
    ## optional ##
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required='run' not in sys.argv)
    parser.add_argument("--reps", type=int, help="The number of replicates to run. [default = 1]", required=False)
    parser.add_argument("--seed", type=str, help="a seed to the random number generator so runs can be replicated", required=False)
    parser.add_argument("--runid", type=str, help="a string to prepend to output files so that multiple runs can be run at the same time without file name clashes", required=False)

    args = parser.parse_args()

    # check --mode
    valid_modes = ["run"]
    if args.mode not in valid_modes:
        print(args.mode, "not a valid mode ( valid modes:",",".join(valid_modes),")", file=sys.stderr)

    #check -r
    args.reference = get_abs_path(args.reference)
    #check -c
    args.consensus = get_abs_path(args.consensus)
    # check -g
    args.locations = get_abs_path(args.locations)
    # check -t
    args.taxonomy = get_abs_path(args.taxonomy)
    # check -j
    args.config = get_abs_path(args.config)
    with open(args.config, "r") as j:
        config = json.load(j, object_pairs_hook = OrderedDict)
    args.config = config

    #check -p
    if args.proc is None:
        args.proc = 1

    #check -o
    if args.out is None:
        args.out = os.path.abspath(".")
    else:
        args.out = os.path.abspath(args.out)

        if not os.path.exists(args.out):
            try:
                os.mkdir(args.out)
            except Exception as e:
                track = traceback.format_exc()
                print(track, file=sys.stderr)
                print("cannot create output directory: ",args.out,"exiting...", file=sys.stderr)
                sys.exit(1)

    # check --reps
    if args.reps is None:
        args.reps = 1
    
    if args.runid is None:
        args.runid = ""

    return args


def writelog(log, msg):
    if log is not None:
        with open(log, "a") as out:
            out.write(msg)

def get_abs_path(in_file, log=None):
    if os.path.isfile(in_file):
        return os.path.abspath(in_file)
    else:
        msg = " ".join(["Cannot find file:",in_file,"exiting....\n"])
        sys.stderr.write(msg)
        writelog(log, msg)
        sys.exit(1)


def get_seqs(fasta):
    seq_dir = []
    records = SeqIO.parse(fasta, "fasta")
    for record in records:
        seq_dir.append([str(record.id), str(record.seq)])
    
    return seq_dir

def fix_fasta_lines(fasta, length):
    lines = []
    fasta_records = SeqIO.parse(fasta,"fasta")
    for record in fasta_records:
        # print(">"+record.id)
        header = ">"+str(record.id)
        lines.append(header)
        seq = str(record.seq)
        x = 0
        while(x+length < len(seq)):
            # print(seq[x:x+length])
            lines.append(seq[x:x+length])
            x += length

        remainder = (len(seq)) - x
        # print(seq[x:x+remainder])
        lines.append(seq[x:x+remainder])
    
    return lines

def add_synthetic_insertion(reference_seqs, consensus_seqs, config, rep, out, run_id="", seed=None):
    if not os.path.exists(out+"/data"):
        os.mkdir(out+"/data")
    
    if seed is not None:
        random.seed(seed+"add_synthetic_insertion"+str(rep))
    else:
        random.seed(str(datetime.now())+"add_synthetic_insertion"+str(rep))

    # get TE family to add
    te_families = list(config["families"].keys())
    family_to_add = te_families[random.randint(0, len(te_families)-1)]
    family_seq = ""
    for seq in consensus_seqs:
        if seq[0] == family_to_add:
            family_seq = seq[1]
    target_file = config['families'][family_to_add]["targets"]
    valid_chroms = []
    with open(target_file,"r") as t:
        for line in t:
            chrom = line.split("\t")[0]
            valid_chroms.append(chrom)
    
    # removed chroms that are not in target file
    tmp = []
    for ref in reference_seqs:
        if ref[0] in valid_chroms:
            tmp.append(ref)
    reference_seqs = tmp

    # get reference contig to modify
    chrom_idx_to_change = random.randint(0,len(reference_seqs)-1)
    chrom_to_change = reference_seqs[chrom_idx_to_change][0]
    seq_to_change = reference_seqs[chrom_idx_to_change][1]


    # get target site to add TE
    targets = []
    with open(target_file,"r") as t:
        for line in t:
            split_line = line.split("\t")
            if split_line[0] == chrom_to_change:
                targets.append(split_line)

    if len(targets) == 0:
        print(chrom_to_change)
        sys.exit("missing chromosome in target files")
    target = targets[random.randint(0, len(targets)-1)]
    target_site = random.randint(int(target[1]), int(target[2])-1)

    duplication_size = config['families'][family_to_add]["TSD"]

    seq_start = seq_to_change[0:target_site+duplication_size]
    seq_end = seq_to_change[target_site:]
    modified_seq = seq_start + family_seq + seq_end

    reference_seqs[chrom_idx_to_change][1] = modified_seq

    with open(out+"/data/"+run_id+str(rep)+".modref.fasta.tmp","w") as outfa:
        for seq in reference_seqs:
            outfa.write(">"+seq[0]+"\n")
            outfa.write(seq[1]+"\n")
    
    fasta_lines = fix_fasta_lines(out+"/data/"+run_id+str(rep)+".modref.fasta.tmp", 80)
    print(len(fasta_lines))
    with open(out+"/data/"+run_id+str(rep)+".modref.fasta","w") as outfa:
        for line in fasta_lines:
            outfa.write(line+"\n")

    os.remove(out+"/data/"+run_id+str(rep)+".modref.fasta.tmp")
    
    with open(out+"/data/"+run_id+str(rep)+".modref.bed","w") as outbed:
        outbed.write(chrom_to_change+"\t"+str(target_site)+"\t"+str(target_site+duplication_size)+"\t"+family_to_add+"\n")

    print(chrom_to_change+"\t"+str(target_site)+"\t"+str(target_site+duplication_size)+"\t"+family_to_add+"\n")

    return out+"/data/"+run_id+str(rep)+".modref.fasta"


if __name__ == "__main__":                
    main()