# create bed file for regions upstream tRNA TSS
# This script is modified based on https://github.com/bergmanlab/preston/blob/master/projects/mcclintock/make_target_sites.py
import sys

gff = sys.argv[1]
upstream_start = int(sys.argv[2]) # 800
upstream_end = int(sys.argv[3]) # 50

with open(gff,"r") as g:
    for line in g:
        split_line = line.split("\t")
        chrom = split_line[0]
        start = int(split_line[3])-1
        end = int(split_line[4])
        strand = split_line[6]

        if strand == "+":
            new_start = start - upstream_start
            new_end = start - upstream_end
        else:
            new_start = end + upstream_end
            new_end = end + upstream_start
        
        if new_start < 0:
            new_start = 0

        if new_end < 0:
            new_end = 0
           
        print(chrom, new_start, new_end, sep="\t")