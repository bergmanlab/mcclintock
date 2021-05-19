# https://github.com/tk2/RetroSeq/wiki/RetroSeq-Tutorial
PARAMS = {
    "read_support_threshold" : 0,
    "breakpoint_confidence_threshold" : 6
}

'''
read_support_threshold                  Number of correctly mapped read pairs spanning breakpoint

breakpoint_confidence_threshold         The FL tag ranges from 1-8 and gives information on the breakpoint with 8 being the most confident calls and 
                                        lower values indicating calls that don’t meet the breakpoint criteria for reasons such as lack of 5’ or 3’ reads
'''