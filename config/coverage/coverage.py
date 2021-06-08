PARAMS = {
    # omits the edges of the TE sequence coverage when calculating average depth across the element
    "omit_edges": True,
    # If OMIT_EDGES = True, use the read length as the length of the edges to omit
    "omit_edges_read_length" : True,
    # IF OMIT_EDGES = True and OMIT_EDGES_READ_LENGTH = False
    # use this value as the length of the edges to omit
    "omit_edges_length" : 300
}
