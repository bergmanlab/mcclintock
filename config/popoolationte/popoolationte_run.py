'''
NAME
    perl identify-te-insertsites.pl - Identifies TE insertion sites (forward
    or reverse insertion) from a sam file
Usage:
     perl identify-te-insertsites.pl --input pe_maped_pool.sam --output te_insertionsites.txt --min-count 3 --narrow-range 100 --min-map-qual 15 --te-hierarchy-file fb-hierarchy.txt --te-hierarchy-level family

    --min-count
        the minimum number of PE-fragments that confirm the insertion of a
        TE of a certain family
    --min-map-qual
        the minimum mapping quality; this will only apply to reads mapping
        to a reference contig.
'''
IDENTIFY_TE_INSERTSITES = {
    "min-count" : 3,
    "min-map-qual" : 15
}




CROSSLINK_TE_SITES = {
    "single-site-shift": 100
}



UPDATE_TEINSERTS_WITH_KNOWNTES = {
    "single-site-shift": 100
}


'''
NAME
    perl estimate-polymorphism.pl - Estimte the insertion frequencies for a
    given set of TE insertions

SYNOPSIS
     perl estimate-polymorphism.pl --sam-file pe_maped_reads.sam --te-insert-file te_insertions.txt --output te-insertion-polymorphism.txt --min-map-qual 15 --te-hierarchy-file te_hierarcy.txt --te-hierarchy-level family

OPTIONS
    --min-map-qual
        the minimum mapping quality

'''

ESTIMATE_POLYMORPHISM = {
    "min-map-qual": 15
}


FILTER = {
    "min-count": 5
}