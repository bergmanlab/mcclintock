
===========
Input files
===========
.. warning::

    Feature names (contig IDs, TE IDs, Family IDs) must not contain any contain any invalid symbols to ensure compatibility with all component methods
    
    INVALID_SYMBOLS :code:`; & ( ) | * ? [ ] ~ { } < ! ^ " , \  $ / + - #`

Required
--------

Reference FASTA (:code:`-r/--reference`)
  * The genome sequence of the reference genome in FASTA format. The reads from the FASTQ file(s) will be mapped to this reference genome to predict TE insertions
  * `Reference FASTA example <https://github.com/bergmanlab/mcclintock/blob/master/test/sacCer2.fasta>`_

Consensus FASTA (:code:`-c/--consensus`)
  * A FASTA file containing a consensus sequence for each TE family.
  * `Consensus FASTA example <https://github.com/bergmanlab/mcclintock/blob/master/test/sac_cer_TE_seqs.fasta>`_

FASTQ 1 (:code:`-1/--first`)
  * Either the Read1 FASTQ file from a paired-end seuqencing run or the FASTQ file from an unpaired sequencing run

Optional
--------

FASTQ 2 (:code:`-2/--second`)
  * The Read2 FASTQ file from a paired-end sequencing run. Not required if using unpaired data.

Locations (:code:`-g/--locations`)
  * A GFF file containing the locations of the reference TEs. 
  * Each annotation should contain an :code:`ID=` attribute that contains a unique identifier for that particular TE. 
  * If this file is provided, the taxonomy file (:code:`-t`) must also be provided.
  * `locations GFF example <https://github.com/bergmanlab/mcclintock/blob/master/test/reference_TE_locations.gff>`_

Taxonomy (:code:`-t/--taxonomy`)
  * A tab-delimited file that maps the unique reference TE ID to the family it belongs to. 
  * This file should contain two columns, the first corresponding to the reference TE identifier which should match the :code:`ID=` attribute from the locations GFF (:code:`-g`).
  * The second column contains the reference TE's family which should match the name of a sequence in the consensus FASTA (:code:`-c`)
  * `Taxonomy example <https://github.com/bergmanlab/mcclintock/blob/master/test/sac_cer_te_families.tsv>`_

Coverage FASTA (:code:`-s/--coverage_fasta`)
  * A FASTA file of TE sequences to be used for the covearge analysis.
  * By default, McClintock estimates the coverage and creates coverage plots of the consensus TE sequences (:code:`-c`). This option allows you to use a custom set of TEs for the coverage estimations and plots.

Augment FASTA (:code:`-a/--augment`)
  * A FASTA file of TE sequences that will be included as extra chromosomes in the reference genome file (:code:`-r`)
  * Some methods leverage the reference TE sequences to find non-reference TEs insertions. The augment FASTA can be used to augment the reference genome with additional TEs that can be used to locate non-reference TE insertions that do not have a representative in the reference genome.