'''
C10: read count for "presence reads"
C11: read count for "absence reads"
C13: genotype for every TE (allele frequency for pooled data, present/absent for haploid, present/absent/heterozygous for diploid) #Note: haploid/diploid caller is under construction, use "pooled" for presence/absence read counts
'''


PARAMETERS = {
    "min_presence_reads" : 3,
    "max_absence_reads" : None,
    "min_presence_fraction" : 0.7 # presence/(absence+presence)
}