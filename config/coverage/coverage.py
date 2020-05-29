
# omits the edges of the TE sequence coverage when calculating average depth across the element
# (start = start + avg_read_length; end = end - avg_read_length)
OMIT_EDGES = True