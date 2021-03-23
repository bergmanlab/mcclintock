rule teflon_preprocessing:
    input:
        te_gff = config['mcc']['unaugmented_locations'],
        taxonomy = config['mcc']['unaugmented_taxonomy'],
        consensus = config['mcc']['consensus'],
        reference_genome = config['mcc']['unaugmented_reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['teflon']

    params:
        out_dir = config['outdir']['teflon']+"unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/teflon/",
        log=config['args']['log_dir']+"teflon.log",
        status_log = config['status']['teflon']

    output:
        ref_bed = config['outdir']['teflon']+"unfiltered/reference_te.bed",
        teflon_taxonomy = config['outdir']['teflon']+"unfiltered/teflon_taxonomy.tsv",
        bam = config['outdir']['teflon']+"unfiltered/teflon.sorted.bam"

    script:
        config['args']['mcc_path']+"/scripts/teflon/teflon_pre.py"


rule teflon_run:
    input:
        consensus = config['mcc']['consensus'],
        reference_genome = config['mcc']['unaugmented_reference'],
        ref_bed = config['outdir']['teflon']+"unfiltered/reference_te.bed",
        teflon_taxonomy = config['outdir']['teflon']+"unfiltered/teflon_taxonomy.tsv",
        bam = config['outdir']['teflon']+"unfiltered/teflon.sorted.bam"

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['teflon']

    params:
        out_dir = config['outdir']['teflon']+"unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/teflon/",
        log=config['args']['log_dir']+"teflon.log",
        config = config['config']['teflon']['files'][0],
        status_log = config['status']['teflon']
    
    output:
        config['outdir']['teflon']+"unfiltered/genotypes/sample.genotypes.txt"
    
    script:
        config['args']['mcc_path']+"/scripts/teflon/teflon_run.py"


rule teflon_post:
    input:
        teflon_out = config['outdir']['teflon']+"unfiltered/genotypes/sample.genotypes.txt",
        ref_bed = config['outdir']['teflon']+"unfiltered/reference_te.bed",
        reference_fasta = config['mcc']['unaugmented_reference']

    threads: 1

    conda: config['envs']['processing']

    params:
        out_dir = config['outdir']['teflon'],
        log=config['args']['log_dir']+"teflon.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        config = config['config']['teflon']['files'][1],
        status_log = config['status']['teflon']
    
    output:
        out = config['out']['teflon']
    
    script:
        config['args']['mcc_path']+"/scripts/teflon/teflon_post.py"
