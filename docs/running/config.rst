
###################
Method config files
###################

Each component method has associated config file(s) that enable the modification of run parameters and filtering parameters. TE detection components have separate run and post processing config files. These files are located in the `config directory <https://github.com/bergmanlab/mcclintock/tree/master/config>`_ of mcclintock :code:`/path/to/mcclintock/config`.

***********************
Example run config file
***********************
.. code:: python

    PARAMS = {
        '-l' : 10,
        '-m' : 0.0,
        '-bm' : 10,
        '-bt' : 7,
        '-f' : 100
    }

Each config file contains a python :code:`dict` named :code:`PARAMS` which contains :code:`key : value` pairs with the :code:`key` being the flag or name of the parameter, and the :code:`value` being the value associated with the parameter. Most of the :code:`values` are set to the default parameters for each tool. The :code:`values` can be modified to better fit the data being used.

********
Coverage
********
:code:`/path/to/mcclintock/config/coverage/coverage.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/coverage/coverage.py>`_

.. code:: python

    PARAMS = {
        "omit_edges": True,
        "omit_edges_read_length" : True,
        "omit_edges_length" : 300
    }

This config file contains the parameters that can be modified for the :code:`coverage` component method.

omit_edges
  * Omits the edges of the TE sequence coverage when calculating average depth across the element.

omit_edges_read_length
  * If :code:`omit_edges: True` and :code:`omit_edges_read_length True`, then the read length will be used as the length of the edges to omit

omit_edges_length
  * If :code:`omit_edges: True` and :code:`omit_edges_read_length False`, the value of :code:`omit_edges_length` will be used as the length of the edges to omit

*************
ngs_te_mapper
*************

run config
==========
:code:`/path/to/mcclintock/config/ngs_te_mapper/ngs_te_mapper_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/ngs_te_mapper/ngs_te_mapper_run.py>`_

.. code:: python

    PARAMS = {
        "tsd=" : 20
    }

This config file contains the parameters that can be modified to influence the running of the :code:`ngs_te_mapper` component method.

tsd=
  * The max length of junction read overlap to consider a target site duplication

postprocessing config
=====================
:code:`/path/to/mcclintock/config/ngs_te_mapper/ngs_te_mapper_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/ngs_te_mapper/ngs_te_mapper_post.py>`_

.. code:: python

    PARAMS = {
        "min_read_support" : 0
    }

This config file contains the parameters that influence the post processing of the :code:`ngs_te_mapper` predictions.

min_read_support
  * the minimum number of reads supporting an insertion required to be considered a valid prediction


**************
ngs_te_mapper2
**************

run config
==========
:code:`/path/to/mcclintock/config/ngs_te_mapper2/ngs_te_mapper2_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/ngs_te_mapper2/ngs_te_mapper2_run.py>`_

.. code:: python

  PARAMS = {
      "--window" : 10,
      "--min_mapq" : 20,
      "--min_af" : 0.1,
      "--tsd_max" : 25,
      "--gap_max" : 5
  }

This config file contains the parameters that can be modified to influence the running of the :code:`ngs_te_mapper2` component method.

--window
  * size of the window to merge for identifying TE clusters

--min_mapq
  * minimum mapping quality of alignment

--min_af
  * minimum allele frequency

--tsd_max
  * maximum target site duplication size

--gap_max
  * maximum gap size

*************
popoolationte
*************

run config
==========
:code:`/path/to/mcclintock/config/popoolationte/popoolationte_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/popoolationte/popoolationte_run.py>`_

.. code:: python

  PARAMS = {
      "identify-te-insertsites.pl" : {
          "--min-count" : 3,
          "--min-map-qual" : 15
      },
      "crosslink-te-sites.pl" : {
          "--single-site-shift": 100
      },
      "update-teinserts-with-knowntes.pl" : {
          "--single-site-shift": 100
      },
      "estimate-polymorphism.pl" : {
          "--min-map-qual": 15
      },
      "filter-teinserts.pl" : {
          "--min-count": 5
      }
  }

identify-te-insertsites.pl
  * Identifies TE insertion sites (forward or reverse insertion) from a sam file
  * :code:`--min-count`: the minimum number of PE-fragments that confirm the insertion of a TE of a certain family
  * :code:`--min-map-qual`: the minimum mapping quality; this will only apply to reads mapping to a reference contig.

crosslink-te-sites.pl
  * Crosslinks forward and reverse insertions and outputs transposable element insertions
  * :code:`--single-site-shift`: the exact position of a TE insertion can only be approximated for TE insertions where only the forward or only the reverse insertion was found. For forward insertions the positon of the TE insertion is calculated as the end of the range. For reverse insertions the position of the TE insertion is calculated as the start of the range

update-teinserts-with-knowntes.pl
  * :code:`--single-site-shift`: the exact position of a TE insertion can only be approximated for TE insertions where only the forward or only the reverse insertion was found. For forward insertions the positon of the TE insertion is calculated as the end of the range. For reverse insertions the position of the TE insertion is calculated as the start of the range

estimate-polymorphism.pl
  * Estimate the insertion frequencies for a given set of TE insertions
  * :code:`--min-map-qual`: the minimum mapping quality

filter-teinserts.pl
  * :code:`--min-count`: the minimum number of PE-fragments that confirm the insertion of a TE of a certain family

postprocessing config
=====================
:code:`/path/to/mcclintock/config/popoolationte/popoolationte_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/popoolationte/popoolationte_post.py>`_

.. code:: python

  PARAMS = {
      "require_both_end_support" : True,
      "percent_read_support_threshold" : 0.1
  }

require_both_end_support
  * requires that final results have support on both ends of prediction

percent_read_support_threshold
  * threshold for the minimum acceptable fraction of the reads supporting the prediction


**************
popoolationte2
**************

run config
==========
:code:`/path/to/mcclintock/config/popoolationte2/popoolationte2_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/popoolationte2/popoolationte2_run.py>`_

.. code:: python

  PARAMS = {
      "ppileup" : {
          "--map-qual": 15,
          "--sr-mindist" : 10000,
          "--id-up-quant": 0.01
      },

      "subsampleppileup" : {
          "run" : False,
          "--target-coverage": 100,
          "--with-replace": False
      },

      "identifySignatures" : {
          "--min-count": 2.0,
          "--signature-window": "median",
          "--min-valley": "median",
          "--chunk-distance": 5
      },

      "updateStrand" : {
          "--map-qual": 15,
          "--max-disagreement": 0.1,
          "--sr-mindist": 10000,
          "--id-up-quant": 0.01
      },

      "pairupSignatures" : {
          "--min-distance": -100,
          "--max-distance": 500,
          "--max-freq-diff": 1.0
      }

  }

ppileup
  * create a physical pileup file from one or multiple bam files
  * :code:`--map-qual`: minimum mapping quality
  * :code:`--sr-mindist`: minimum inner distance for structural rearrangements
  * :code:`--id-up-quant`: paired end fragments with an insert size in the upper quantile will be ignored

subsampleppileup
  * subsample a ppileup file to uniform coverage
  * :code:`run`: The subsampleppileup step is optional. Set this option to :code:`True` if you wish to perform this step
  * :code:`--target-coverage`: the target coverage of the output file
  * :code:`--with-replace`: use sampling with replacement instead of without replacement

identifySignatures
  * identify signatures of TE insertions
  * :code:`--min-count`: the minimum count of a TE insertion
  * :code:`--signature-window`: the window size of the signatures of TE insertions. options:  :code:`median`, :code:`fixNNNN`, :code:`minimumSampleMedian`, :code:`maximumSampleMedian`
  * :code:`min-valley`: the minimum size of the valley between two consecutive signatures of the same family. options:  :code:`median`, :code:`fixNNNN`, :code:`minimumSampleMedian`, :code:`maximumSampleMedian`
  * :code:`--chunk-distance`: minimum distance between chromosomal chunks in multiples of :code:`--min-valley`

updateStrand
  * estimate the strand of TEs for signatures of TE insertions
  * :code:`--map-qual`: minimum mapping quality
  * :code:`--max-disagreement`: the maximum disagreement for the strand of the TE insertion in fraction of reads
  * :code:`--sr-mindist`: minimum inner distance for structural rearrangements
  * :code:`--id-up-quant`: paired end fragments with an insert size in the upper quantile will be ignored

pairupSignatures
  * pairs up signatures of TE insertions and yields TE insertions
  * :code:`--min-distance`: the minimum distance between signatures
  * :code:`--max-distance`: the maximum distance between signatures
  * :code:`--max-freq-diff`: the maximum frequency difference between signatures


postprocessing config
=====================
:code:`/path/to/mcclintock/config/popoolationte2/popoolationte2_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/popoolationte2/popoolationte2_post.py>`_

.. code:: python

  PARAMS = {
      "require_both_end_support" : True,
      "frequency_threshold" : 0.1
  }

require_both_end_support
  * require that the TE prediction have support from both junctions

frequency_threshold
  * threshold for the minimum acceptable (average physical coverage supporting the insertion of the given TE) / (average physical coverage)


********
relocate
********

run config
==========
:code:`/path/to/mcclintock/config/relocate/relocate_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/relocate/relocate_run.py>`_

.. code:: python

  PARAMS = {
      '-l' : 10,
      '-m' : 0.0,
      '-bm' : 10,
      '-bt' : 7,
      '-f' : 100
  }

-l
  * len cutoff for the TE trimmed reads to be aligned

-m
  * mismatch allowance for alignment to TE

-bm
  * blat minScore value, used by blat in the comparison of reads to TE sequence

-bt
  * blat tileSize value, used by blat in the comparison of reads to TE sequence

-f
  * length of the sequence flanking the found insertion to be returned. This sequence is taken from the reference genome


postprocessing config
=====================
:code:`/path/to/mcclintock/config/relocate/relocate_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/relocate/relocate_post.py>`_

.. code:: python

  PARAMS = {
      "ref_left_threshold" : 0,
      "ref_right_threshold" : 0,
      "nonref_left_threshold" : 0,
      "nonref_right_threshold" : 0
  }

ref_left_threshold
  * minimum number of left flanking reads for a reference prediction.

ref_right_threshold
  * minimum number of right flanking reads for a reference prediction.

nonref_left_threshold
  * minimum number of left flanking reads for a non-reference prediction.

nonref_right_threshold
  * minimum number of right flanking reads for a non-reference prediction.


*********
relocate2
*********

run config
==========
:code:`/path/to/mcclintock/config/relocate2/relocate2_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/relocate2/relocate2_run.py>`_

.. code:: python

  PARAMS = {
      '--aligner' : "blat",
      '--len_cut_match' : 10,
      '--len_cut_trim' : 10,
      '--mismatch' : 2,
      '--mismatch_junction' : 2
  }

--aligner
  * aligner used to map reads to repeat elements

--len_cut_match
  * length cutoff threshold for match between reads and repeat elements. Large value will lead to less sensitive but more accuracy

--len_cut_trim
  * length cutoff threshold for trimed reads after trimming repeat sequence from reads. Large value will lead to less sensitive but more accuracy

--mismatch
  * Number of mismatches allowed for matches between reads and repeat elements

--mismatch_junction
  * Number of mismatches allowed for matches between junction reads and repeat elements


postprocessing config
=====================
:code:`/path/to/mcclintock/config/relocate2/relocate2_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/relocate2/relocate2_post.py>`_

.. code:: python

  PARAMS = {
      "ref_left_support_threshold" : 0,
      "ref_right_support_threshold" : 0,
      "ref_left_junction_threshold" : 0,
      "ref_right_junction_threshold" : 0,

      "nonref_left_support_threshold" : 0,
      "nonref_right_support_threshold" : 0,
      "nonref_left_junction_threshold" : 0,
      "nonref_right_junction_threshold" : 0
  }

ref_left_support_threshold
  * Minimum number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on left side/downstream of a reference prediction

ref_right_support_threshold
  * Minimum number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on right side/downstream of a reference prediction

ref_left_junction_threshold
  * Minimum number of reads covering the junction of TE insertion on left side/upstream of a reference prediction

ref_right_junction_threshold
  * Minimum number of reads covering the junction of TE insertion on right side/downstream of a reference prediction

nonref_left_support_threshold
  * Minimum number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on left side/downstream of a non-reference prediction

nonref_right_support_threshold
  * Minimum number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on right side/downstream of a non-reference prediction

nonref_left_junction_threshold
  * Minimum number of reads covering the junction of TE insertion on left side/upstream of a non-reference prediction

nonref_right_junction_threshold
  * Minimum number of reads covering the junction of TE insertion on right side/downstream of a non-reference prediction


********
retroseq
********

run config
==========
:code:`/path/to/mcclintock/config/retroseq/retroseq_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/retroseq/retroseq_run.py>`_

.. code:: python

  PARAMS = {
      "-depth" : 200,
      "-reads" : 10,
      "-q": 20
  }

-depth
  * Max average depth of a region to be considered for calling

-reads
  * Minimum number of reads required to make a call

-q
  * Minimum mapping quality for a read mate that anchors the insertion call


postprocessing config
=====================
:code:`/path/to/mcclintock/config/retroseq/retroseq_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/retroseq/retroseq_post.py>`_

.. code:: python

  PARAMS = {
      "read_support_threshold" : 0,
      "breakpoint_confidence_threshold" : 6
  }

read_support_threshold
  * Minimum number of correctly mapped read pairs spanning breakpoint for predictions

breakpoint_confidence_threshold
  * Minimum FL tag for predictions. The FL tag ranges from 1-8 and gives information on the breakpoint with 8 being the most confident calls and lower values indicating calls that don’t meet the breakpoint criteria for reasons such as lack of 5’ or 3’ reads


*******
tebreak
*******

run config
==========
:code:`/path/to/mcclintock/config/tebreak/tebreak_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/tebreak/tebreak_run.py>`_

.. code:: python

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

--minMWP
  * minimum Mann-Whitney P-value for split qualities

--min_minclip
  * min. shortest clipped bases per cluster

--min_maxclip
  * min. longest clipped bases per cluster

--min_sr_per_break
  * minimum split reads per breakend

--min_consensus_score
  * quality of consensus alignment

--min_chr_len
  * minimum chromosome length to consider in discordant site search

--max_ins_reads
  * maximum number of reads to use per insertion call

--min_split_reads
  * minimum total split reads per insertion call

--min_prox_mapq
  * minimum map quality for proximal subread

--max_N_consensus
  * exclude breakend seqs with > this number of N bases

--max_disc_fetch
  * maximum number of discordant reads to fetch per insertion site per BAM

--min_disc_reads
  * if using -d/--disco_target, minimum number of discordant reads to trigger a call

--sr_density
  * maximum split read density in chunk

--min_ins_match
  * (output) minumum match to insertion library

--min_ref_match
  * (output) minimum match to reference genome

--min_cons_len
  * (output) min total consensus length

--keep_all_tmp_bams
  * leave ALL temporary BAMs

--skip_final_filter
  * do not apply final filters or fix for orientation

--debug
  * run in debug mode

postprocessing config
=====================
:code:`/path/to/mcclintock/config/tebreak/tebreak_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/tebreak/tebreak_post.py>`_

.. code:: python

  MIN_5P_ELT_MATCH = 0.0
  MIN_3P_ELT_MATCH = 0.0
  MIN_5P_GENOME_MATCH = 0.0
  MIN_3P_GENOME_MATCH = 0.0
  MIN_SPLIT_READS_5P = 0
  MIN_SPLIT_READS_3P = 0
  MIN_REMAPPED_DISCORDANT = 0
  MIN_REMAP_DISC_FRACTION = 0.0
  MIN_REMAPPED_SPLITREADS = 0
  MIN_REMAP_SPLIT_FRACTION = 0.0

MIN_5P_ELT_MATCH
  * Minimum fraction of bases matched to reference for inserted sequence on insertion seqment of 5' supporting contig

MIN_3P_ELT_MATCH
  * Minimum fraction of bases matched to reference for inserted sequence on insertion seqment of 3' supporting contig

MIN_5P_GENOME_MATCH
  * Minimum fraction of bases matched to reference genome on genomic segment of 5' supporting contig

MIN_3P_GENOME_MATCH
  * Minimum fraction of bases matched to reference genome on genomic segment of 3' supporting contig

MIN_SPLIT_READS_5P
  * Minimum number of split reads supporting 5' end of the insertion

MIN_SPLIT_READS_3P
  * Minimum number of split reads supporting 3' end of the insertion

MIN_REMAPPED_DISCORDANT
  * Minimum number of discordant read ends re-mappable to insertion reference sequence

MIN_REMAP_DISC_FRACTION
  * Minimum proportion of remapped discordant reads mapping to the reference insertion sequence

MIN_REMAPPED_SPLITREADS
  * Minimum number of split reads re-mappable to insertion reference sequence

MIN_REMAP_SPLIT_FRACTION
  * Minimum proportion of remapped split reads mapping to the reference insertion sequence

******
teflon
******

run config
==========
:code:`/path/to/mcclintock/config/teflon/teflon_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/teflon/teflon_run.py>`_

.. code:: python

  PARAMS = {
      "-q": 20,
      "-sd" : None,
      "-cov" : None,
      "-n1" : 1,
      "-n2" : 1,
      "-lt" : 1,
      "-ht" : None
  }

-q
  * Map quality threshold. Mapped reads with map qualities lower than this number will be discarded

-sd
  * Insert size standard deviation. Used to manually override the insert size StdDev identified by samtools stat (check this number in the generated stats.txt file to ensure it seems more or less correct based on knowledge of sequencing library!)
  * Set to :code:`None` if you want this value to be calculated by TEFLoN

-cov
  * Coverage override. Used to manually override the coverage estimate if you get the error: :code:`Warning: coverage could not be estimated`

-n1
  * TEs must be supported by >= n reads in at least one sample

-n2
  * TEs must be supported by >= n reads summed across all samples

-lt
  * sites genotyped as -9 if adjusted read counts lower than this threshold

-ht
  * sites genotyped as -9 if adjusted read counts higher than this threshold
  * Set to :code:`None` if you want this value to be calculated by TEFLoN


postprocessing config
=====================
:code:`/path/to/mcclintock/config/teflon/teflon_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/teflon/teflon_post.py>`_

.. code:: python

  PARAMS = {
      "min_presence_reads" : 3,
      "max_absence_reads" : None,
      "min_presence_fraction" : 0.1,
      "require_tsd" : True,
      "require_both_breakpoints" : True

  }

min_presence_reads
  * Minimum number of reads supporting the presence of the insertion

max_absence_reads
  * Maximum number of reads that support the absence of the insertion

min_presence_fraction
  * The minimum fraction of reads supporting the presence of the insertion. presence/(absence+presence)

require_tsd
  * If :code:`True`, non-ref predictions must have a target site duplication

require_both_breakpoints
  * If :code:`True`, non-ref predictions must have both breakpoints predicted


*********
te-locate
*********

run config
==========
:code:`/path/to/mcclintock/config/telocate/telocate_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/telocate/telocate_run.py>`_

.. code:: python

  PARAMS = {
      "max_mem" : 4,
      "min_distance" : 5,
      "min_support_reads" : 3,
      "min_support_individuals" : 1
  }

max_mem
  * Max memory to use (GB)

min_distance
  * resolution for the loci, if a supporting read-pair is found in a distance up to this value it is counted for the same event. Multiplied by median insert size.

min_support_reads
  * only events supported by this number of reads in all accessions are kept

min_support_individuals
  * only events supported in this number of individuals are kept.

postprocessing config
=====================
:code:`/path/to/mcclintock/config/telocate/telocate_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/telocate/telocate_post.py>`_

.. code:: python

  PARAMS = {
      "read_pair_support_threshold" : 0
  }

read_pair_support_threshold
  * Minimum number of read pairs supporting a prediction

****
temp
****

run config
==========
:code:`/path/to/mcclintock/config/temp/temp_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/temp/temp_run.py>`_

.. code:: python

  PARAMS = {
      '-x' : 30,
      '-m' : 1
  }

-x
  * The minimum score difference between the best hit and the second best hit for considering a read as uniquely mapped. For BWA mem.

-m
  * Number of mismatch allowed when mapping to TE concensus sequences.


postprocessing config
=====================
:code:`/path/to/mcclintock/config/temp/temp_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/temp/temp_post.py>`_

.. code:: python

  PARAMS = {
      "acceptable_insertion_support_classes" : ["1p1"],
      "frequency_threshold" : 0.1
  }

acceptable_insertion_support_classes
  * valid options: :code:`1p1`, :code:`2p`, :code:`singleton`
  * The class of the insertions that are acceptable for final predictions. 
  * :code:`1p1` means that the detected insertion is supported by reads at both sides. 
  * :code:`2p` means the detected insertion is supported by more than 1 read at only 1 side. 
  * :code:`singleton` means the detected insertion is supported by only 1 read at 1 side.

frequency_threshold
  * minimum frequency of the inserted transposon


*****
temp2
*****

run config
==========
:code:`/path/to/mcclintock/config/temp2/temp2_run.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/temp2/temp2_run.py>`_

.. code:: python

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

insertion
  * :code:`-M` : Percentage of mismatch allowed when anchor to genome.
  * :code:`-m` : Percentage of mismatch allowed when mapping to TEs.
  * :code:`-U` : The ratio between the second best alignment and the best alignment to judge if a read is uniquely mapped.
  * :code:`-N` : window size (+-n) for filtering insertions overlapping reference insertions.
  * :code:`-T` : Set this parameter to :code:`True` to allow truncated de novo insertions.
  * :code:`-L` : Set this parameter to :code:`True` to use a looser criteria to filter reference annotated copy overlapped insertions.
  * :code:`-S` : Set this parameter to :code:`True` to skip insertion length checking; Default is to remove those insertions that are not full length of shorter than 500bp.

absence
  * :code:`-x` : The minimum score difference between the best hit and the second best hit for considering a read as uniquely mapped. For BWA MEM.


postprocessing config
=====================
:code:`/path/to/mcclintock/config/temp2/temp2_post.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/temp2/temp2_post.py>`_

.. code:: python

  PARAMS = {
      "acceptable_insertion_support_classes" : ["1p1"],
      "frequency_threshold" : 0.1
  }

acceptable_insertion_support_classes
  * valid options: :code:`1p1`, :code:`2p`, :code:`singleton`
  * The class of the insertions that are acceptable for final predictions. 
  * :code:`1p1` means that the detected insertion is supported by reads at both sides. 
  * :code:`2p` means the detected insertion is supported by more than 1 read at only 1 side. 
  * :code:`singleton` means the detected insertion is supported by only 1 read at 1 side.

frequency_threshold
  * minimum frequency of the inserted transposon. It generally means what fraction of sequenced genome present this insertion.


**********
trimgalore
**********

:code:`/path/to/mcclintock/config/trimgalore/trimgalore.py` : `GitHub link <https://github.com/bergmanlab/mcclintock/blob/master/config/trimgalore/trimgalore.py>`_

.. code:: python

  PARAMS = {
      "single_end": {
          "--fastqc": True,
      },

      "paired_end": {
          "--fastqc": True,
          "--paired": True
      }
  }

single_end
  * :code:`--fastqc`: Set to true to run fastqc when using single end data

paired_end
  * :code:`--fastqc`: Set to true to run fastqc when using paired end data
  * :code:`--paired`: Run trimgalore using length trimming of quality/adapter/RRBS trimmed reads for paired-end data