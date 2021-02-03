
##############################
# TEMP2 INSERTION PARAMETERS #
##############################

# -M mismatch%    Percentage of mismatch allowed when anchor to genome. Default is 2.
GENOME_MISMATCH_PCT = 2

# -m mismatch%    Percentage of mismatch allowed when mapping to TEs. Default is 5.
TE_MISMATCH_PCT = 5

# -U ratio    The ratio between the second best alignment and the best alignment to judge if a read is uniquely mapped. Default is 0.8.
RATIO = 0.8

# -N reference_filter_window      window sizea (+-n) for filtering insertions overlapping reference insertions. Default is 300.
FILTER_WINDOW = 300

# -T              Set this parameter to allow truncated de novo insertions; For default, only full-length de novo insertions are allowed.
TRUNCATED = False

# -L              Set this parameter to use a looser criteria to filter reference annotated copy overlapped insertions; Default not allowed.
LOOSE_FILTER = False

# -S              Set this parameter to skip insertion length checking; Default is to remove those insertions that are not full length of shorter than 500bp.
SKIP_INS_LEN_CHECK = False


############################
# TEMP2 ABSENCE PARAMETERS #
############################

# -x     The minimum score difference between the best hit and the second best hit for considering a read as uniquely mapped. For BWA MEM.
UNIQ_MAP_SCORE = 0