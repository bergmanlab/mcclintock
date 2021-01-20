
================
Running examples
================

Here are some examples of how to run McClintock for common use cases.

Run McClintock with test data
-----------------------------
Some test data is provided in the :code:`test/` directory of the `McClintock repository <https://github.com/bergmanlab/mcclintock>`_. FastQ files must be downloaded using the :code:`test/download_test_data.py` script.

.. code:: bash

    python test/download_test_data.py

The test data provided is: 

* a UCSC sacCer2 yeast reference genome ( :code:`sacCer2.fasta` )
* consensus sequences of the TY families in sacCer ( :code:`sac_cer_TE_seqs.fasta` )
* an annotation of TEs in the yeast reference genome from `Carr, Bensasson and Bergman (2012) <http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0050978>`_ ( :code:`reference_TE_locations.gff` )
* A taxonomy file linking reference TE ids from :code:`reference_TE_locations.gff` to the TE family they belong to. ( :code:`sac_cer_te_families.tsv` )
* a pair of fastq files from NCBI Sequence Read Archive ( :code:`SRR800842_1.fastq.gz` and :code:`SRR800842_2.fastq.gz` )

To run mcclintock on the test data, use the following command:

.. code:: bash

    python3 mcclintock.py \
        -r test/sacCer2.fasta \
        -c test/sac_cer_TE_seqs.fasta \
        -g test/reference_TE_locations.gff \
        -t test/sac_cer_te_families.tsv \
        -1 test/SRR800842_1.fastq.gz \
        -2 test/SRR800842_2.fastq.gz \
        -p 4 \
        -o /path/to/output/directory

change :code:`/path/to/output/directory` to a real path where you desire the McClintock output to be created. You can also increase :code:`-p 4` to a higher nuber if you have more CPU threads available.

Run McClintock with specific component methods
----------------------------------------------
By default, McClintock runs all component methods with the data provided. If you only want to run a specific component method, you can use the :code:`-m` flag to specify which method to run

.. code:: bash

    python3 mcclintock.py \
        -r test/sacCer2.fasta \
        -c test/sac_cer_TE_seqs.fasta \
        -g test/reference_TE_locations.gff \
        -t test/sac_cer_te_families.tsv \
        -1 test/SRR800842_1.fastq.gz \
        -2 test/SRR800842_2.fastq.gz \
        -p 4 \
        -m temp \
        -o /path/to/output/directory

You can also specify multiple methods to run by writing a comma-separated list of the methods after the :code:`-m` flag

.. code:: bash

    python3 mcclintock.py \
        -r test/sacCer2.fasta \
        -c test/sac_cer_TE_seqs.fasta \
        -g test/reference_TE_locations.gff \
        -t test/sac_cer_te_families.tsv \
        -1 test/SRR800842_1.fastq.gz \
        -2 test/SRR800842_2.fastq.gz \
        -p 4 \
        -m temp,ngs_te_mapper,retroseq \
        -o /path/to/output/directory

Run McClintock with multiple samples using same reference genome
----------------------------------------------------------------
When running McClintock on multiple samples that use the same reference genome and consensus TEs, it is advised to pre-generate the TE locations GFF and a TE Taxonomy TSV. Otherwise, these files will be redundantly created by mcclintock for each sample. If you lack a TE locations GFF and a TE Taxonomy TSV, you can run McClintock with the :code:`--make_annotations` flag to produce these files in advance.

.. code:: bash

    python3 mcclintock.py \
        -r test/sacCer2.fasta \
        -c test/sac_cer_TE_seqs.fasta \
        -p 4 \
        -o <output> \
        --make_annotations

With the :code:`--make_annotations` flag, McClintock will produce the reference TE locations GFF and taxonomy file using RepeatMasker, then exit the run.

* Reference TE locations GFF: :code:`<output>/<reference_name>/reference_te_locations/unaugmented_inrefTEs.gff`
* TE Taxonomy TSV: :code:`<output>/<reference_name>/te_taxonomy/unaugmented_taxonomy.tsv`

You can then use the :code:`--resume` flag for future runs with the same reference genome and output directory without having to redundantly generate them for each run.

.. code:: bash

    python3 mcclintock.py \
        -r test/sacCer2.fasta \
        -c test/sac_cer_TE_seqs.fasta \
        -1 /path/to/sample1_1.fastq.gz \
        -2 /path/to/sample1_2.fastq.gz \
        -p 4 \
        -o <output> \
        --resume

    python3 mcclintock.py \
        -r test/sacCer2.fasta \
        -c test/sac_cer_TE_seqs.fasta \
        -1 /path/to/sample2_1.fastq.gz \
        -2 /path/to/sample2_2.fastq.gz \
        -p 4 \
        -o <output> \
        --resume

    ## etc ##

Individual samples can be run in a serial manner as shown in the example above, or run in parallel, such as through separate jobs on a HPC cluster.