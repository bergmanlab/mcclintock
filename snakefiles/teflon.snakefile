rule teflon_preprocessing:
    input:
        te_gff = config['mcc']['locations'],
        taxonomy = config['mcc']['taxonomy'],
        consensus = config['mcc']['consensus'],
        reference_genome = config['mcc']['reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['teflon']

    params:
        out_dir = config['args']['out']+"results/teflon/unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/teflon/",
        log=config['args']['log_dir']+"teflon.log"

    output:
        ref_bed = config['args']['out']+"results/teflon/unfiltered/reference_te.bed",
        teflon_taxonomy = config['args']['out']+"results/teflon/unfiltered/teflon_taxonomy.tsv",
        bam = config['args']['out']+"results/teflon/unfiltered/teflon.sorted.bam"

    script:
        config['args']['mcc_path']+"/scripts/teflon/teflon_pre.py"


rule teflon_run:
    input:
        consensus = config['mcc']['consensus'],
        reference_genome = config['mcc']['reference'],
        ref_bed = config['args']['out']+"results/teflon/unfiltered/reference_te.bed",
        teflon_taxonomy = config['args']['out']+"results/teflon/unfiltered/teflon_taxonomy.tsv",
        bam = config['args']['out']+"results/teflon/unfiltered/teflon.sorted.bam"

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['teflon']

    params:
        out_dir = config['args']['out']+"results/teflon/unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/teflon/",
        log=config['args']['log_dir']+"teflon.log"
    
    output:
        config['args']['out']+"results/teflon/unfiltered/genotypes/sample.genotypes.txt"
    
    script:
        config['args']['mcc_path']+"/scripts/teflon/teflon_run.py"


rule teflon_post:
    input:
        teflon_out = config['args']['out']+"results/teflon/unfiltered/genotypes/sample.genotypes.txt",
        ref_bed = config['args']['out']+"results/teflon/unfiltered/reference_te.bed"

    threads: 1

    conda: config['envs']['teflon']

    params:
        out_dir = config['args']['out']+"results/teflon/",
        log=config['args']['log_dir']+"teflon.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']
    
    output:
        out = config['out']['teflon']
    
    script:
        config['args']['mcc_path']+"/scripts/teflon/teflon_post.py"
