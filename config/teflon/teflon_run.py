
# When params: sd, cov, or ht are set to None, they are determined by TEFLoN, Can be set to a specific INT value 
PARAMS = {
    "-q": 20,
    "-sd" : None,
    "-cov" : None,
    "-n1" : 1,
    "-n2" : 1,
    "-lt" : 1,
    "-ht" : None
}

'''
usage: python usr/local/teflon.v0.4.py <required> [optional]
    -q <map quality threshold> #NOTE: Mapped reads with map qualities lower than this number will be discarded
    -sd [insert size standard deviation] #NOTE: Used to manually override the insert size StdDev identified by samtools stat (check this number in the generated stats.txt file to ensure it seems more or less correct based on knowledge of sequencing library!)
    -cov [coverage override] #Note: Used to manually override the coverage estimate if you get the error: "Warning: coverage could not be estimated"

usage: python usr/local/teflon_collapse.py <required> [optional]
    -n1 <TEs must be supported by >= n reads in at least one sample>
    -n2 <TEs must be supported by >= n reads summed across all samples>
    -q <map quality threshold> #NOTE: Mapped reads with map qualities lower than this number will be discarded
    -cov [coverage override] #NOTE: Used to manually override the coverage estimate if you get the error: "Warning: coverage could not be estimated"

usage: python /usr/local/teflon_count.py <required> [optional]
    -q <map quality threshold>

usage: python usr/local/teflon_genotype.py <required> [optional]
    -lt [sites genotyped as -9 if adjusted read counts lower than this threshold, default=1]
    -ht [sites genotyped as -9 if adjusted read counts higher than this threshold, default=mean_coverage + 2*STDEV]
'''
