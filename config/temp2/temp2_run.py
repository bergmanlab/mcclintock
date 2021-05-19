PARAMS ={
    "insertion" : {
        "-M" : 2,
        "-m" : 5,
        "-U" : 0.8,
        "-N" : 300,
        "-T" : False,
        "-L" : False,
        "-S" : False
    },

    "absence" : {
        "-x" : 0
    }
}

##############################
# TEMP2 INSERTION PARAMETERS #
##############################

# -M mismatch%    Percentage of mismatch allowed when anchor to genome. Default is 2.
# -m mismatch%    Percentage of mismatch allowed when mapping to TEs. Default is 5.
# -U ratio    The ratio between the second best alignment and the best alignment to judge if a read is uniquely mapped. Default is 0.8.
# -N reference_filter_window      window sizea (+-n) for filtering insertions overlapping reference insertions. Default is 300.
# -T              Set this parameter to allow truncated de novo insertions; For default, only full-length de novo insertions are allowed.
# -L              Set this parameter to use a looser criteria to filter reference annotated copy overlapped insertions; Default not allowed.
# -S              Set this parameter to skip insertion length checking; Default is to remove those insertions that are not full length of shorter than 500bp.


############################
# TEMP2 ABSENCE PARAMETERS #
############################

# -x     The minimum score difference between the best hit and the second best hit for considering a read as uniquely mapped. For BWA MEM.