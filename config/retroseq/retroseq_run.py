'''
[-depth         Max average depth of a region to be considered for calling. Default is 200.]
[-reads         It is the minimum number of reads required to make a call. Default is 5. Parameter is optional.]
[-q             Minimum mapping quality for a read mate that anchors the insertion call. Default is 30. Parameter is optional.]
'''

PARAMETERS = {
    "depth" : 200,
    "reads" : 10,
    "q": 20
}