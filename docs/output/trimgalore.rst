
==========
TrimGalore
==========

The TrimGalore module runs `Trim Galore <https://github.com/FelixKrueger/TrimGalore>`_ on the input FastQ files for QC and to trim adaptors prior to running the component methods. The output of the trimgalore module can be found in :code:`<output>/<sample>/results/trimgalore/`.

:code:`<fastq>_trimming_report.txt`

* Information on parameters used and statistics related to adapter trimming with cutadapt. Provides an overview of sequences removed via the adapter trimming process.

:code:`<fastq>_fastqc.html`

* FastQC report of the trimmed fastq files. Provides information on the results of steps performed by FastQC to assess the quality of the trimmed reads.

:code:`<fastq>_fastqc.zip`

* FastQC graph images and plain-text summary statistics compressed into a single :code:`.zip` file