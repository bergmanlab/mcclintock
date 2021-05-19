
PARAMS = {
    "min_presence_reads" : 3,
    "max_absence_reads" : None,
    "min_presence_fraction" : 0.1, # presence/(absence+presence),
    "require_tsd" : True, # non-ref predictions must have a target site duplication
    "require_both_breakpoints" : True # non-ref predictions must have both breakpoints predicted

}

'''
C10: read count for "presence reads"
C11: read count for "absence reads"
C13: genotype for every TE (allele frequency for pooled data, present/absent for haploid, present/absent/heterozygous for diploid) #Note: haploid/diploid caller is under construction, use "pooled" for presence/absence read counts
'''

