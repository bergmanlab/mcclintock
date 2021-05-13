rule tebreak_run:
    input:
        consensus_fasta = config['mcc']['consensus'],
        bam = config['mcc']['bam'],
        ref_fasta = config['mcc']['reference'],
        rm_out = config['mcc']['repeatmasker_out']

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['tebreak']

    params:
        script_dir = config['args']['mcc_path']+"/install/tools/tebreak/tebreak/",
        out_dir = config['outdir']['tebreak']+"unfiltered/",
        ref_name=config['args']['ref_name'],
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"tebreak.log",
        config = config['config']['tebreak']['files'][0],
        status_log = config['status']['tebreak']
    
    output:
        config['outdir']['tebreak']+"unfiltered/"+config['args']['sample_name']+".sorted.tebreak.table.txt"
    
    script:
        config['args']['mcc_path']+"/scripts/tebreak/tebreak_run.py"

rule tebreak_post:
    input:
        tebreak_out = config['outdir']['tebreak']+"unfiltered/"+config['args']['sample_name']+".sorted.tebreak.table.txt",
        ref_fasta = config['mcc']['reference']

    threads: 1

    conda: config['envs']['processing']

    params:
        out_dir = config['outdir']['tebreak'],
        ref_name=config['args']['ref_name'],
        sample_name=config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        config = config['config']['tebreak']['files'][1],
        status_log = config['status']['tebreak']
    
    output:
        config['out']['tebreak']
    
    script:
        config['args']['mcc_path']+"/scripts/tebreak/tebreak_post.py"