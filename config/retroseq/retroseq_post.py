# https://github.com/tk2/RetroSeq/wiki/RetroSeq-Tutorial
'''
The calls are annotated with information on number of supporting reads (GQ tag)
'''
READ_SUPPORT_THRESHOLD = 0

'''
The FL tag ranges from 1-8 and gives information on the breakpoint with 8 being the most confident calls and 
lower values indicating calls that don’t meet the breakpoint criteria for reasons such as lack of 5’ or 3’ reads
'''
BREAKPOINT_CONFIDENCE_THRESHOLD = 6