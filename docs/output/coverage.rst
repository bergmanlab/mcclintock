
========
Coverage
========

The coverage module estimates copy number based on normalized coverage and creates coverage plots for each TE in the fasta provided by :code:`-c/--consensus` or :code:`-s/--coverage_fasta` if provided. The coverage profiles can be viewed as images or as tables of the raw/unique coverage at each position. The coverage module output files can be found in the :code:`<output>/<sample>/results/trimgalore/` output directory.


:code:`plots/*.png`

* Coverage plots showing the normalized read coverage across each TE either from the consensus fasta (:code:`-c`) or the coverage fasta (:code:`-s`) if provided. Coverage of uniquely mapping reads (MAPQ > 0) is in dark gray, while coverage of all reads (MAPQ >= 0) is in light gray. Raw coverage at each postion in a TE is normalized to the average mapping depth at unique regions of the hard-masked reference genome. The average normalized coverage is shown as a black line, and is estimated from the central region of each TE omitting regions at the 5' and 3' ends equal to the average read length to prevent biases due to mapping at TE edges.

:code:`te-depth-files/*.allQ.cov`

* Raw read coverage at each position in a TE sequence. (Output of :code:`samtools depth`)

:code:`te-depth-files/*.highQ.cov`

* coverage of mapped reads with MAPQ > 0 at each position, omitting multi-mapped reads.