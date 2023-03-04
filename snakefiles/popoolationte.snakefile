rule make_popoolationte_annotations:
    input:
        te_gff = config['mcc']['unaugmented_locations'],
        taxonomy = config['mcc']['unaugmented_taxonomy'],
        consensus = config['mcc']['consensus']

    params:
        mcc_out = config['args']['out'],
        run_id = config['args']['run_id'],
        log = config['args']['log_dir']+"processing.log",
        chromosomes = config['args']['chromosomes'],
        augment = config['args']['augment_fasta']

    threads: 1

    conda: config['envs']['processing']

    output:
        taxonomy = config['mcc']['popoolationTE_taxonomy'],
        te_gff = config['mcc']['popoolationTE_gff']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_popoolationte_annotations.py"


rule popoolationTE_ref_fasta:
    input:
        config['mcc']['masked_fasta'],
        config['mcc']['consensus'],
        config['mcc']['ref_te_fasta']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['processing']

    output:
        config['mcc']['popoolationTE_ref_fasta']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/popoolationTE_ref_fasta.py"

rule popoolationTE_preprocessing:
    input:
        ref_fasta = config['mcc']['popoolationTE_ref_fasta'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['popoolationte']

    params:
        out_dir = config['outdir']['popoolationte']+"unfiltered/",
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"popoolationTE.log",
        script_dir = config['args']['mcc_path']+"/install/tools/popoolationte/",
        status_log = config['status']['popoolationte']

    output:
        config['outdir']['popoolationte']+"unfiltered/reads1.fastq",
        config['outdir']['popoolationte']+"unfiltered/reads2.fastq",
        config['outdir']['popoolationte']+"unfiltered/combined.sorted.sam"

    script:
        config['args']['mcc_path']+"/scripts/popoolationte/popoolationte_pre.py"

rule popoolationTE_run:
    input:
        ref_fasta = config['mcc']['popoolationTE_ref_fasta'],
        taxonomy = config['mcc']['popoolationTE_taxonomy'],
        te_gff = config['mcc']['popoolationTE_gff'],
        fq1 = config['outdir']['popoolationte']+"unfiltered/reads1.fastq",
        fq2 = config['outdir']['popoolationte']+"unfiltered/reads2.fastq",
        sam = config['outdir']['popoolationte']+"unfiltered/combined.sorted.sam"
    
    threads: 1

    conda: config['envs']['popoolationte']

    params:
        out_dir = config['outdir']['popoolationte']+"unfiltered/",
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"popoolationTE.log",
        script_dir = config['args']['mcc_path']+"/install/tools/popoolationte/",
        config = config['config']['popoolationte']['files'][0],
        status_log = config['status']['popoolationte']

    output:
        config['outdir']['popoolationte']+"unfiltered/te-poly-filtered.txt"

    script:
        config['args']['mcc_path']+"/scripts/popoolationte/popoolationte_run.py"

rule popoolationTE_post:
    input:
        popoolationte_out = config['outdir']['popoolationte']+"unfiltered/te-poly-filtered.txt",
        ref = config['mcc']['unaugmented_reference']
    
    threads: 1

    conda: config['envs']['processing']

    params:
        out_dir = config['outdir']['popoolationte'],
        sample_name=config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        log = config['args']['log_dir']+"popoolationTE.log",
        config = config['config']['popoolationte']['files'][1],
        status_log = config['status']['popoolationte'],
        vcf = config['args']['vcf']

    output:
        config['out']['popoolationte']

    script:
        config['args']['mcc_path']+"/scripts/popoolationte/popoolationte_post.py"