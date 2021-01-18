rule ngs_te_mapper2_run:
    input:
        config = config['args']['mcc_path']+"/config/ngs_te_mapper2/ngs_te_mapper2_run.py",
        consensus_fasta = config['mcc']['consensus'],
        reference_fasta = config['mcc']['reference'],
        fastq1 = config['mcc']['fq1'],
        fastq2 = config['mcc']['fq2']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/ngs_te_mapper2/unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/ngs_te_mapper2/sourceCode/",
        log = config['args']['log_dir']+"ngs_te_mapper2.log",
        sample_name = config['args']['sample_name']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['ngs_te_mapper2']

    output:
        config['args']['out']+"results/ngs_te_mapper2/unfiltered/"+config['args']['sample_name']+".nonref.bed",
        config['args']['out']+"results/ngs_te_mapper2/unfiltered/"+config['args']['sample_name']+".ref.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/ngs_te_mapper2/ngs_te_mapper2_run.py"


rule ngs_te_mapper2_post:
    input:
        config = config['args']['mcc_path']+"/config/ngs_te_mapper2/ngs_te_mapper2_post.py",
        ref_bed = config['args']['out']+"results/ngs_te_mapper2/unfiltered/"+config['args']['sample_name']+".ref.bed",
        nonref_bed = config['args']['out']+"results/ngs_te_mapper2/unfiltered/"+config['args']['sample_name']+".nonref.bed",
        reference_fasta = config['mcc']['reference']

    params:
        out_dir = config['args']['out']+"/results/ngs_te_mapper2/",
        log = config['args']['log_dir']+"ngs_te_mapper2.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']
    
    threads: 1

    conda: config['envs']['mcc_processing']

    output:
        config['out']['ngs_te_mapper2']
    
    script:
        config['args']['mcc_path']+"/scripts/ngs_te_mapper2/ngs_te_mapper2_post.py"