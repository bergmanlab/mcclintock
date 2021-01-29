rule ngs_te_mapper_run:
    input:
        config = config['args']['mcc_path']+"/config/ngs_te_mapper/ngs_te_mapper_run.py",
        consensus_fasta = config['mcc']['consensus'],
        reference_fasta = config['mcc']['reference'],
        fastq1 = config['mcc']['fq1'],
        fastq2 = config['mcc']['fq2']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/ngs_te_mapper/unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/ngs_te_mapper/sourceCode/",
        log = config['args']['log_dir']+"ngs_te_mapper.log",
        sample_name = config['args']['sample_name']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['ngs_te_mapper']

    output:
        config['args']['out']+"results/ngs_te_mapper/unfiltered/bed_tsd/"+config['args']['sample_name']+"_insertions.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/ngs_te_mapper/ngs_te_mapper_run.py"

rule ngs_te_mapper_post:
    input:
        config = config['args']['mcc_path']+"/config/ngs_te_mapper/ngs_te_mapper_post.py",
        raw_bed = config['args']['out']+"results/ngs_te_mapper/unfiltered/bed_tsd/"+config['args']['sample_name']+"_insertions.bed",
        reference_fasta = config['mcc']['reference']

    params:
        out_dir = config['args']['out']+"/results/ngs_te_mapper/",
        log = config['args']['log_dir']+"ngs_te_mapper.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']
    
    threads: 1

    conda: config['envs']['processing']

    output:
        config['out']['ngs_te_mapper']
    
    script:
        config['args']['mcc_path']+"/scripts/ngs_te_mapper/ngs_te_mapper_post.py"