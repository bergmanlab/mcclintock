rule relocaTE2_run:
    input:
        config['args']['mcc_path']+"/config/relocate2/relocate2_run.py",
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
        out_dir = config['args']['out']+"/results/relocaTE2/unfiltered/",
        log = config['args']['log_dir']+"relocaTE2.log",
        sample_name = config['args']['sample_name']
        
    
    output:
        config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_nonref_insert.gff",
        config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_ref_insert.gff"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE2/relocate2_run.py"



rule relocaTE2_post:
    input:
        config['args']['mcc_path']+"/config/relocate2/relocate2_post.py",
        rm_out = config['mcc']['repeatmasker_out'],
        nonref_gff = config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_nonref_insert.gff",
        ref_gff = config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_ref_insert.gff"

    threads: 1

    conda: config['envs']['relocate2']

    params:
        out_dir = config['args']['out']+"/results/relocaTE2/",
        log = config['args']['log_dir']+"relocaTE2.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']
    
    output:
        config['args']['out']+"/results/relocaTE2/"+config['args']['sample_name']+"_relocate2_redundant.bed",
        config['args']['out']+"/results/relocaTE2/"+config['args']['sample_name']+"_relocate2_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE2/relocate2_post.py"