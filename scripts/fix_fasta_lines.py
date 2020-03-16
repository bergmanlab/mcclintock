#!/usr/bin/env python3

import sys
from Bio import SeqIO


def main():
    fasta = sys.argv[1]
    length = int(sys.argv[2])

    fasta_records = SeqIO.parse(fasta,"fasta")
    for record in fasta_records:
        print(">"+record.id)
        seq = str(record.seq)
        x = 0
        while(x+length < len(seq)):
            print(seq[x:x+length])
            x += length

        remainder = (len(seq)) - x
        print(seq[x:x+remainder])


if __name__ == "__main__":                
    main()