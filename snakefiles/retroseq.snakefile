rule retroseq_run:
    input:
        consensus_fasta = config['mcc']['consensus'],
        bam = config['mcc']['bam'],
        ref_fasta = config['mcc']['reference'],
        ref_te_bed = config['mcc']['ref_tes_bed'],
        taxonomy = config['mcc']['taxonomy']

    threads: 1

    conda: config['envs']['retroseq']

    params:
        script_dir = config['args']['mcc_path']+"/install/tools/retroseq/bin/",
        out_dir = config['args']['out']+"/results/retroseq/unfiltered/",
        ref_name=config['args']['ref_name'],
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"retroseq.log"
    
    output:
        config['args']['out']+"results/retroseq/unfiltered/"+config['args']['sample_name']+".call"
    
    script:
        config['args']['mcc_path']+"/scripts/retroseq/retroseq_run.py"

rule retroseq_post:
    input:
        retroseq_out = config['args']['out']+"results/retroseq/unfiltered/"+config['args']['sample_name']+".call",
        reference_fasta = config['mcc']['reference']

    threads: 1

    conda: config['envs']['mcc_processing']

    params:
        out_dir = config['args']['out']+"/results/retroseq/",
        ref_name=config['args']['ref_name'],
        sample_name=config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']
    
    output:
        config['out']['retroseq']
    
    script:
        config['args']['mcc_path']+"/scripts/retroseq/retroseq_post.py"