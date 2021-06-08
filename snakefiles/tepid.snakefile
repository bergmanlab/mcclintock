rule tepid_run:
    input:
        ref_fasta = config['mcc']['reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        te_gff = config['mcc']['locations'],
        te_taxonomy = config['mcc']['taxonomy'],
        median_insert_size = config['mcc']['median_insert_size']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['tepid']

    params:
        raw_fq2 = config['in']['fq2'],
        ref_name = config['args']['ref_name'],
        out_dir = config['args']['out']+"/results/tepid/unfiltered/",
        log = config['args']['log_dir']+"tepid.log",
        config = config['config']['tepid']['files'][0],
        status_log = config['status']['tepid']
    
    output:
        config['args']['out']+"results/tepid/unfiltered/insertions_"+config['args']['ref_name']+".bed",
        config['args']['out']+"results/tepid/unfiltered/deletions_"+config['args']['ref_name']+".bed",
        insertions_support = config['args']['out']+"results/tepid/unfiltered/insertion_reads_"+config['args']['ref_name']+".txt",
        deletions_support = config['args']['out']+"results/tepid/unfiltered/deletion_reads_"+config['args']['ref_name']+".txt"

    script:
        config['args']['mcc_path']+"/scripts/tepid/tepid_run.py"

rule tepid_post:
    input:
        insertions_bed = config['args']['out']+"results/tepid/unfiltered/insertions_"+config['args']['ref_name']+".bed",
        deletions_bed = config['args']['out']+"results/tepid/unfiltered/deletions_"+config['args']['ref_name']+".bed",
        insertions_support = config['args']['out']+"results/tepid/unfiltered/insertion_reads_"+config['args']['ref_name']+".txt",
        deletions_support = config['args']['out']+"results/tepid/unfiltered/deletion_reads_"+config['args']['ref_name']+".txt",
        te_gff = config['mcc']['locations'],
        te_taxonomy = config['mcc']['taxonomy'],
        reference_fasta = config['mcc']['reference']

    threads: 1

    conda: config['envs']['tepid']

    params:
        sample_name = config['args']['sample_name'],
        out_dir = config['args']['out']+"/results/tepid/",
        chromosomes = config['args']['chromosomes'],
        config = config['config']['tepid']['files'][1],
        status_log = config['status']['tepid']

    output:
        config['out']['tepid']
    
    script:
        config['args']['mcc_path']+"/scripts/tepid/tepid_post.py"