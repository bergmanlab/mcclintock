PARAMS = {
    "--minMWP": "0.01",
    "--min_minclip" : "3",
    "--min_maxclip" : "10",
    "--min_sr_per_break" : "1",
    "--min_consensus_score" : "0.9",
    "--min_chr_len" : "0",
    "--max_ins_reads" : "1000",
    "--min_split_reads" : "4",
    "--min_prox_mapq" : "10",
    "--max_N_consensus" : "4",
    "--max_disc_fetch" : "50",
    "--min_disc_reads" : "4",
    "--sr_density" : "2.0",
    "--min_ins_match" : "0.90",
    "--min_ref_match" : "0.98",
    "--min_cons_len" : "250",
    "--keep_all_tmp_bams" : False,
    "--skip_final_filter" : False,
    "--debug": False
}


#   --minMWP MINMWP                            minimum Mann-Whitney P-value for split qualities (default = 0.01)
#   --min_minclip MIN_MINCLIP                  min. shortest clipped bases per cluster (default = 3)
#   --min_maxclip MIN_MAXCLIP                  min. longest clipped bases per cluster (default = 10)
#   --min_sr_per_break MIN_SR_PER_BREAK        minimum split reads per breakend (default = 1)
#   --min_consensus_score MIN_CONSENSUS_SCORE  quality of consensus alignment (default = 0.9)
#   --min_chr_len MIN_CHR_LEN                  minimum chromosome length to consider in discordant site search (default = 0)
#   --max_ins_reads MAX_INS_READS              maximum number of reads to use per insertion call (default = 1000)
#   --min_split_reads MIN_SPLIT_READS          minimum total split reads per insertion call (default = 4)
#   --min_prox_mapq MIN_PROX_MAPQ              minimum map quality for proximal subread (default = 10)
#   --max_N_consensus MAX_N_CONSENSUS          exclude breakend seqs with > this number of N bases (default = 4)
#   --max_disc_fetch MAX_DISC_FETCH            maximum number of discordant reads to fetch per insertion site per BAM (default = 50; 0 = disable fetch)
#   --min_disc_reads MIN_DISC_READS            if using -d/--disco_target, minimum number of discordant reads to trigger a call (default = 4)
#   --sr_density SR_DENSITY                    maximum split read density in chunk (default = 2.0)
#   --min_ins_match MIN_INS_MATCH              (output) minumum match to insertion library (default 0.90)
#   --min_ref_match MIN_REF_MATCH              (output) minimum match to reference genome (default 0.98)
#   --min_cons_len MIN_CONS_LEN                (output) min total consensus length (default=250)
#   --keep_all_tmp_bams                        leave ALL temporary BAMs (warning: lots of files!)
#   --skip_final_filter                        do not apply final filters or fix for orientation
#   --debug