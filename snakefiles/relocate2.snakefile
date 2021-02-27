rule relocaTE2_run:
    input:
        reference = config['mcc']['reference'],
        te_seqs = config['mcc']['consensus'],
        rm_out = config['mcc']['repeatmasker_out'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        median_insert_size = config['mcc']['median_insert_size']

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['relocate2']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['outdir']['relocate2']+"unfiltered/",
        log = config['args']['log_dir']+"relocaTE2.log",
        sample_name = config['args']['sample_name'],
        config = config['config']['relocate2']['files'][0]
        
    
    output:
        config['outdir']['relocate2']+"unfiltered/repeat/results/ALL.all_nonref_insert.gff",
        config['outdir']['relocate2']+"unfiltered/repeat/results/ALL.all_ref_insert.gff"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE2/relocate2_run.py"



rule relocaTE2_post:
    input:
        rm_out = config['mcc']['repeatmasker_out'],
        nonref_gff = config['outdir']['relocate2']+"unfiltered/repeat/results/ALL.all_nonref_insert.gff",
        ref_gff = config['outdir']['relocate2']+"unfiltered/repeat/results/ALL.all_ref_insert.gff",
        reference_fasta = config['mcc']['reference']

    threads: 1

    conda: config['envs']['processing']

    params:
        out_dir = config['outdir']['relocate2'],
        log = config['args']['log_dir']+"relocaTE2.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        config = config['config']['relocate2']['files'][1]
    
    output:
        config['out']['relocate2']
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE2/relocate2_post.py"