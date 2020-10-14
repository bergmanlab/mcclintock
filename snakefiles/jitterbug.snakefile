rule jitterbug_run:
    input:
        reference_tes = config['mcc']['locations'],
        bam = config['mcc']['bam']

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['jitterbug']

    params:
        out_dir = config['args']['out']+"results/jitterbug/unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/jitterbug/",
        sample_name = config['args']['sample_name'],
        log=config['args']['log_dir']+"jitterbug.log"
    
    output:
        out = config['args']['out']+"results/jitterbug/unfiltered/"+config['args']['sample_name']+".TE_insertions_paired_clusters.filtered.gff3"
    
    script:
        config['args']['mcc_path']+"/scripts/jitterbug/jitterbug_run.py"

rule jitterbug_post:
    input:
        jitterbug_out = config['args']['out']+"results/jitterbug/unfiltered/"+config['args']['sample_name']+".TE_insertions_paired_clusters.filtered.gff3",
        taxonomy = config['mcc']['taxonomy']

    threads: 1

    conda: config['envs']['jitterbug']

    params:
        out_dir = config['args']['out']+"results/jitterbug/",
        log=config['args']['log_dir']+"jitterbug.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']
    
    output:
        out = config['out']['jitterbug']
    
    script:
        config['args']['mcc_path']+"/scripts/jitterbug/jitterbug_post.py"