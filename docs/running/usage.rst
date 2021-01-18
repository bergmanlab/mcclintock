
=====
Usage
=====
The main McClintock script can be run in three primary modes: `run`, `install`, and `make_annotations`

Run usage
---------
* McClintock, when run without the :code:`--install` and :code:`--make_annotations` flags, will run the complete TE annotation prediction pipeline.

.. code:: text

    usage: mcclintock.py [-h] -r REFERENCE -c CONSENSUS -1 FIRST [-2 SECOND]
                    [-p PROC] [-o OUT] [-m METHODS] [-g LOCATIONS] [-t TAXONOMY]
                    [-s COVERAGE_FASTA] [-T] [-a AUGMENT]
                    [--sample_name SAMPLE_NAME] [--resume] [--debug]
                    [--slow] [-k KEEP_INTERMEDIATE]

    Meta-pipeline to identify transposable element insertions using next
    generation sequencing data

    required arguments:
    -r REFERENCE, --reference REFERENCE
                            A reference genome sequence in fasta format
    -c CONSENSUS, --consensus CONSENSUS
                            The consensus sequences of the TEs for the species in
                            fasta format
    -1 FIRST, --first FIRST
                            The path of the first fastq file from paired end read
                            sequencing or the fastq file from single read
                            sequencing

    optional arguments:
    -h, --help            show this help message and exit
    -2 SECOND, --second SECOND
                            The path of the second fastq file from a paired end
                            read sequencing
    -p PROC, --proc PROC  The number of processors to use for parallel stages of
                            the pipeline [default = 1]
    -o OUT, --out OUT     An output folder for the run. [default = '.']
    -m METHODS, --methods METHODS
                            A comma-delimited list containing the software you
                            want the pipeline to use for analysis. e.g. '-m
                            relocate,TEMP,ngs_te_mapper' will launch only those
                            three methods
    -g LOCATIONS, --locations LOCATIONS
                            The locations of known TEs in the reference genome in
                            GFF 3 format. This must include a unique ID attribute
                            for every entry
    -t TAXONOMY, --taxonomy TAXONOMY
                            A tab delimited file with one entry per ID in the GFF
                            file and two columns: the first containing the ID and
                            the second containing the TE family it belongs to. The
                            family should correspond to the names of the sequences
                            in the consensus fasta file
    -s COVERAGE_FASTA, --coverage_fasta COVERAGE_FASTA
                            A fasta file that will be used for TE-based coverage
                            analysis, if not supplied then the consensus sequences
                            of the TEs will be used for the analysis
    -T, --comments        If this option is specified then fastq comments (e.g.
                            barcode) will be incorporated to SAM output. Warning:
                            do not use this option if the input fastq files do not
                            have comments
    -a AUGMENT, --augment AUGMENT
                            A fasta file of TE sequences that will be included as
                            extra chromosomes in the reference file (useful if the
                            organism is known to have TEs that are not present in
                            the reference strain)
    --sample_name SAMPLE_NAME
                            The sample name to use for output files [default:
                            fastq1 name]
    --resume              This option will attempt to use existing intermediate
                            files from a previous McClintock run
    --debug               This option will allow snakemake to print progress to
                            stdout
    --slow                This option runs without attempting to optimize thread
                            usage to run rules concurrently. Each multithread rule
                            will use the max processors designated by -p/--proc
    -k KEEP_INTERMEDIATE, --keep_intermediate KEEP_INTERMEDIATE
                            This option determines which intermediate files are
                            preserved after McClintock completes [default:
                            general][options: minimal, general, methods,
                            <list,of,methods>, all]

Install usage
-------------
   McClintock executed with the :code:`--install` flag will run in `install` mode. This will install the component methods and their conda environments. See :doc:`../getting_started/install` for more details.

.. code:: text

    usage: mcclintock.py --install [-h] [-m METHODS] [--resume] [--debug]

    required arguments:
    --install             This option will install the dependencies of
                            mcclintock

    optional arguments:
    -h, --help            show this help message and exit
    -m METHODS, --methods METHODS
                          A comma-delimited list containing the software you
                            want to install. e.g. '-m
                            relocate,TEMP,ngs_te_mapper' will install those
                            three methods and create their conda environments
    --resume              This option will not delete the existing installation,
                            but will install the methods that have not already been installed
    --debug               This option will allow snakemake to print progress to
                            stdout

Make annotations usage
----------------------
   McClintock executed with the :code:`--make_annotations` flag will create the reference TE annotations. This step is part of the main mcclintock pipeline, but this option is useful for creating reference files in advance for future sample runs that will use the same reference genome. See :doc:`examples` for more details on how to use this option.

.. code:: text

    usage: mcclintock.py --make_annotations [-h] -r REFERENCE -c CONSENSUS
                  [-p PROC] [-o OUT] [--resume] [--debug]

    required arguments:
    --make_annotations    This option will only run the pipeline up to the
                            creation of the repeat annotations
    -r REFERENCE, --reference REFERENCE
                            A reference genome sequence in fasta format
    -c CONSENSUS, --consensus CONSENSUS
                            The consensus sequences of the TEs for the species in
                            fasta format

    optional arguments:
    -h, --help            show this help message and exit
    -g LOCATIONS, --locations LOCATIONS
                            The locations of known TEs in the reference genome in
                            GFF 3 format. This must include a unique ID attribute
                            for every entry
    -t TAXONOMY, --taxonomy TAXONOMY
                            A tab delimited file with one entry per ID in the GFF
                            file and two columns: the first containing the ID and
                            the second containing the TE family it belongs to. The
                            family should correspond to the names of the sequences
                            in the consensus fasta file
    -p PROC, --proc PROC  The number of processors to use for parallel stages of
                            the pipeline [default = 1]
    -o OUT, --out OUT     An output folder for the run. [default = '.']
    -a AUGMENT, --augment AUGMENT
                            A fasta file of TE sequences that will be included as
                            extra chromosomes in the reference file (useful if the
                            organism is known to have TEs that are not present in
                            the reference strain)
    --resume              This option will attempt to use existing intermediate
                            files from a previous McClintock run
    --debug               This option will allow snakemake to print progress to
                            stdout