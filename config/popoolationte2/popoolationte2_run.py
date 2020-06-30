
#ppileup
# create a physical pileup file from one or multiple bam files

# == Main parameters ==
# --map-qual            minimum mapping quality; default=15

# == Parameters for fine tuning ==
# --sr-mindist          minimum inner distance for structural rearrangements;
#                       default=10000
# --id-up-quant         paired end fragments with an insert size in the upper
#                       quantile will be ignored [fraction]; default=0.01

ppileup = {
    "map-qual": 15,
    "sr-mindist" : 10000,
    "id-up-quant": 0.01
}


# subsampleppileup 
# NOTE: subsampling of ppileup does not run by default, change "run" : True to turn on this step
# subsample a ppileup file to uniform coverage

# == Main parameters ==
# --target-coverage     the target coverage of the output file [int]; mandatory

# == Parameters for fine tuning ==
# --with-replace        use sampling with replacement instead of without replacement

subsampleppileup = {
    "run" : False,
    "target-coverage": 100,
    "with-replace": False
}


# identifySignatures
# identify signatures of TE insertions

# == Main parameters ==
# --min-count           the minimum count of a TE insertion; default=2.0

# == Parameters for fine tuning ==
# --signature-window    the window size of the signatures of TE insertions;
#                       [median|fixNNNN|minimumSampleMedian|maximumSampleMedian]
#                       ; default=median
# --min-valley          the minimum size of the valley between two consecutive
#                       signatures of the same family ; [median|fixNNNN|minimumSampleMedian|maximumSampleMedian]
#                       ; default=median
# --chunk-distance      minimum distance between chromosomal chunks in multiples
#                       of --min-valley [int]; default=5

identifySignatures = {
    "min-count": 2.0,
    "signature-window": "median",
    "min-valley": "median",
    "chunk-distance": 5
}


# updateStrand
# estimate the strand of TEs for signatures of TE insertions

# == Main parameters ==
# --map-qual            minimum mapping quality; default=15
# --max-disagreement    the maximum disagreement for the strand of the TE insertion
#                       in fraction of reads; mandatory

# == Parameters for fine tuning ==
# --sr-mindist          minimum inner distance for structural rearrangements;
#                       default=10000
# --id-up-quant         paired end fragments with an insert size in the upper
#                       quantile will be ignored [fraction]; default=0.01

updateStrand = {
    "map-qual": 15,
    "max-disagreement": 0.1,
    "sr-mindist": 10000,
    "id-up-quant": 0.01
}


# pairupSignatures
# pairs up signatures of TE insertions and yields TE insertions

# == Main parameters ==
# --min-distance        the minimum distance between signatures; default=-100
# --max-distance        the maximum distance between signatures; default=500

# == Parameters for fine tuning ==
# --max-freq-diff       the maximum frequency difference between signatures;
#                       default=1.0

pairupSignatures = {
    "min-distance": -100,
    "max-distance": 500,
    "max-freq-diff": 1.0
}
