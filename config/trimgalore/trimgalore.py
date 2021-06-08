# options used for adapter trimming with the trimgalore module
PARAMS = {
    "single_end": {
        "--fastqc": True,
    },

    "paired_end": {
        "--fastqc": True,
        "--paired": True
    }
}

