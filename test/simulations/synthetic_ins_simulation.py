import argparse
import os
import sys
import traceback
import subprocess
from multiprocessing import Process, Pool
from Bio import SeqIO


class ModifiedChrom:
    def __init__(self):
        self.name = ""
        self.seq = ""
        self.insert = ""
        self.instert_start = -1
        self.insert_end = -1

class Feature:
    def __init__(self):
        self.chrom = ""
        self.start = -1
        self.end = -1
        self.strand = "?"

def main():
    args = parse_args()
    
    trnas = get_trnas(args.trna)
    chroms_with_inserts = make_inserts(args.reference, args.consensus, trnas)
    if not os.path.exists(args.out+"/data"):
        os.mkdir(args.out+"/data")
    fastas = make_fastas(args.reference, chroms_with_inserts, args.out+"/data", start=args.start, end=args.end)






def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock', description="Meta-pipeline to identify transposable element insertions using next generation sequencing data")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required=True)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required=True)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required=True)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required=True)
    parser.add_argument("--trna", type=str, help="A gff indicating the locations of tRNAs in the reference genome", required=True)
    

    ## optional ##
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required=False)
    parser.add_argument("-s", "--start", type=int, help="the start of the range of seeds to run (default=1)", required=False)
    parser.add_argument("-e", "--end", type=int, help="the end of the range of seeds to run (default=299)", required=False)

    args = parser.parse_args()


    #check -r
    args.reference = get_abs_path(args.reference)
    #check -c
    args.consensus = get_abs_path(args.consensus)
    # check -g
    args.locations = get_abs_path(args.locations)
    # check -t
    args.taxonomy = get_abs_path(args.taxonomy)
    # check --trna
    args.trna = get_abs_path(args.trna)

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

    if args.start is None:
        args.start = 1
    if args.end is None:
        args.end = 299

    if args.start >= args.end:
        print("-s/--start must be lower than -e/--end")
        sys.exit(1)
    
    return args


def get_trnas(trna_gff):
    trnas = []
    with open(trna_gff, "r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                trna = Feature()
                trna.chrom = split_line[0]
                trna.start = int(split_line[3])
                trna.end = int(split_line[4])
                trna.strand = split_line[6]

                trnas.append(trna)
    
    return trnas

def make_inserts(reference, consensus, trnas):
    ref_records = SeqIO.parse(reference, "fasta")
    ref_chroms = {}
    for record in ref_records:
        ref_chroms[str(record.id)] = str(record.seq)

    consensus_records = SeqIO.parse(consensus, "fasta")
    
    consensus_contigs = {}
    for record in consensus_records:
        consensus_contigs[str(record.id)] = str(record.seq)
    
    
    te_types = ["TY1", "TY2", "TY3", "TY4"]
    te_idx = 0
    synthetic_insertions = []
    for trna in trnas:
        if te_idx >= len(te_types):
            te_idx = 0

        insertion = ModifiedChrom()
        seq = ref_chroms[trna.chrom]
        if trna.strand == "+":
            if te_types[te_idx] == "TY3":
                seq_start = seq[0:trna.start-12]
                seq_end = seq[trna.start-17:len(seq)]
            else:
                seq_start = seq[0:trna.start-195]
                seq_end = seq[trna.start-200:len(seq)]
        
        else:
            if te_types[te_idx] == "TY3":
                seq_start = seq[0:trna.start+17]
                seq_end = seq[trna.start+12:len(seq)]
            else:
                seq_start = seq[0:trna.start+200]
                seq_end = seq[trna.start+195:len(seq)]
                # print(seq_start[-5:], seq_end[0:5], len(seq_start+seq_end) - len(seq))
        
        
        insertion.name = trna.chrom
        insertion.seq = seq_start+consensus_contigs[te_types[te_idx]]+seq_end
        insertion.start = len(seq_start)
        insertion.end = len(seq_start+consensus_contigs[te_types[te_idx]])
        synthetic_insertions.append(insertion)
        te_idx += 1

    return synthetic_insertions


def get_abs_path(in_file, log=None):
    if os.path.isfile(in_file):
        return os.path.abspath(in_file)
    else:
        msg = " ".join(["Cannot find file:",in_file,"exiting....\n"])
        sys.stderr.write(msg)
        writelog(log, msg)
        sys.exit(1)


def writelog(log, msg):
    if log is not None:
        with open(log, "a") as out:
            out.write(msg)

def make_fastas(reference, chroms_with_inserts, out, start=1, end=299):
    fastas = []
    for x in range(start-1, end):
        out_fasta = out+"/insertion"+str(x+1)+"_genome.fasta"
        with open(out_fasta, "w") as outfa:
            ref_records = SeqIO.parse(reference, "fasta")
            for record in ref_records:
                if str(record.id) == chroms_with_inserts[x].name:
                    print("replacing chrom with chrom with synthetic insert", x+1)
                    seq = chroms_with_inserts[x].seq
                else:
                    seq = str(record.seq)
                
                outfa.write(">"+str(record.id)+"\n")
                outfa.write(seq+"\n")
        
        fastas.append(out_fasta)
    
    return fastas


if __name__ == "__main__":                
    main()