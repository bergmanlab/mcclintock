
'''
python scripts/relocaTE2.py --help
usage: relocaTE2.py [-h] [-b BAM] [-t TE_FASTA] [-d FQ_DIR] [-g GENOME_FASTA]
                   [-r REFERENCE_INS] [-o OUTDIR] [-s SIZE] [-c CPU]
                   [-1 MATE_1_ID] [-2 MATE_2_ID] [-u UNPAIRED_ID]
                   [--sample SAMPLE] [--aligner ALIGNER]
                   [--len_cut_match LEN_CUT_MATCH]
                   [--len_cut_trim LEN_CUT_TRIM] [--mismatch MISMATCH]
                   [--mismatch_junction MISMATCH_JUNCTION] [--step STEP]
                   [--run] [--split] [-v VERBOSE]

  --aligner ALIGNER     aligner used to map reads to repeat elements,
                        default=blat
  --len_cut_match LEN_CUT_MATCH
                        length cutoff threshold for match between reads and
                        repeat elements. Large value will lead to less
                        sensitive but more accuracy, default = 10
  --len_cut_trim LEN_CUT_TRIM
                        length cutoff threshold for trimed reads after
                        trimming repeat sequence from reads. Large value will
                        lead to less sensitive but more accuracy, default = 10
  --mismatch MISMATCH   Number of mismatches allowed for matches between reads
                        and repeat elements, default = 2
  --mismatch_junction MISMATCH_JUNCTION
                        Number of mismatches allowed for matches between
                        junction reads and repeat elements, default = 2

'''



RELOCATE2 = {
    'aligner' : "blat",
    'len_cut_match' : 10,
    'len_cut_trim' : 10,
    'mismatch' : 2,
    'mismatch_junction' : 2
}