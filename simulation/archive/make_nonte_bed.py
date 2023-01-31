import sys
import argparse
import os
import subprocess
from Bio import SeqIO

class Transposon:
    def __init__(self, chrom, chrom_length, start, end, buffer):
        self.chrom = chrom
        self.chrom_length = chrom_length
        self.start = start
        self.end  = end
        self.buffered_start = max(start-buffer,1)
        self.buffered_end = min(end+buffer, self.chrom_length)


def main():
    args = parse_args()
    chrom_lengths = get_chrom_lengths(args.reference)
    transposons = read_transposons(args.gff, chrom_lengths, buffer=args.buffer)
    te_bed = make_te_bed(transposons, args.out)
    chrom_bed = make_chrom_bed(chrom_lengths, args.out, telomere=args.telomere, omit=args.chroms_to_omit)
    make_nonte_bed(te_bed, chrom_bed, args.out)

    os.remove(te_bed)
    os.remove(chrom_bed)


def parse_args():
    parser = argparse.ArgumentParser(prog='make_nonte_bed', description="Makes a bed file of the non-TE regions in the genome. Useful as targets for McClintock simulation")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required=True)
    parser.add_argument("-g", "--gff", type=str, help="The repeatmasker GFF of repetitive regions in genome", required=True)

    ## optional ##
    parser.add_argument("-o", "--out", type=str, help="File name of the output bed [default=nonTE.bed]", required=False)
    parser.add_argument("-b", "--buffer", type=int, help="buffer region around TEs to omit from final nonTE bed intervals [default=200]", required=False)
    parser.add_argument("-t", "--telomere", type=int, help="How much of the telomeric regions to omit from nonTE bed [default=20000]", required=False)
    parser.add_argument("-c", "--chroms_to_omit", type=str, help="Comma separated list of chromosomes to omit from the nonTE bed [default=None]", required=False)
    
    args = parser.parse_args()

    if args.out is None:
        args.out = "nonTE.bed"
    
    if args.buffer is None:
        args.buffer = 200
    
    if args.telomere is None:
        args.telomere = 20000
    
    if args.chroms_to_omit is None:
        args.chroms_to_omit = []

    return args

def get_chrom_lengths(fasta):
    chrom_lengths = {}
    for record in SeqIO.parse(fasta, "fasta"):
        seqid = str(record.id)
        seq_length = len(str(record.seq))
        chrom_lengths[seqid] = seq_length
    
    return chrom_lengths

def read_transposons(rm_gff, chrom_lengths, buffer=200):
    transposons = []

    with open(rm_gff,"r") as gff:
        for line in gff:
            line = line.strip()
            if "#" not in line:
                split_line = line.split("\t")
                chrom = split_line[0]
                chrom_length = chrom_lengths[chrom]
                start = int(split_line[3])
                end = int(split_line[4])
                tpsn = Transposon(chrom, chrom_length, start, end, buffer)
                transposons.append(tpsn)
    
    return transposons

def make_te_bed(transposons, out):
    out_bed = out+".inverse.tmp.bed"
    with open(out_bed,"w") as o:
        for transposon in transposons:
            o.write("\t".join([transposon.chrom, str(transposon.buffered_start-1),str(transposon.buffered_end),".",".","."])+"\n")
    
    return out_bed

def make_chrom_bed(chrom_lengths, out, telomere=0, omit=[]):
    out_bed = out+".chrom.tmp.bed"
    with open(out_bed, "w") as o:
        for chrom, length in chrom_lengths.items():
            if chrom not in omit:
                start = telomere
                end = length - telomere
                if end < start:
                    sys.exit("Error: telomere length is too long for chromosome:"+chrom+"..exiting..\n")
                
                o.write( "\t".join([chrom, str(start), str(end), ".",".","."])+"\n")
    
    return out_bed

def make_nonte_bed(te_bed, chrom_bed, out):
    out_bed = open(out,"w")
    subprocess.call(["bedtools", "subtract", "-a", chrom_bed, "-b", te_bed], stdout=out_bed)
    out_bed.close()


if __name__ == "__main__":
    main()