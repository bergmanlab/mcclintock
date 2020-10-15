'''
  -q MINMAPQ, --minMAPQ MINMAPQ
                        minimum read mapping quality to be considered
  -c MIN_CLUSTER_SIZE, --min_cluster_size MIN_CLUSTER_SIZE
                        min number of both fwd and rev reads to predict an
                        insertion
'''

RUN = {
    "MINMAPQ" : 15,
    "MIN_CLUSTER_SIZE": 2
}


# Set values to None if you want the value to be determined by jitterbug
FILTER = {
  "CLUSTER_SIZE" : None,
  "SPAN" : None,
  "INT_SIZE" : [0,None],
  "SOFTCLIPPED": [2,None],
  "PICK_CONSISTENT" : None
}

# FILTER = {
#   "CLUSTER_SIZE" : [2, 1687],
#   "SPAN" : [2, 302],
#   "INT_SIZE" : [101, 522],
#   "SOFTCLIPPED": [2, 1687],
#   "PICK_CONSISTENT" : [0, -1]
# }