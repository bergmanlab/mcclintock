rule relocaTE_consensus:
    input:
        config['mcc']['consensus']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['processing']

    output:
        config['mcc']['relocaTE_consensus']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/relocaTE_consensus.py"
    
rule relocaTE_ref_gff:
    input:
        config['mcc']['locations'],
        config['mcc']['taxonomy']

    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['processing']

    output:
        config['mcc']['relocaTE_ref_TEs']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/relocaTE_ref_gff.py" 

rule relocaTE_run:
    input:
        consensus_fasta = config['mcc']['relocaTE_consensus'],
        te_gff = config['mcc']['relocaTE_ref_TEs'],
        reference_fasta = config['mcc']['reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']

    threads: 1

    conda: config['envs']['relocate']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['outdir']['relocate']+"unfiltered/",
        log = config['args']['log_dir']+"relocaTE.log",
        script_dir = config['args']['mcc_path']+"/install/tools/relocate/scripts/",
        sample_name = config['args']['sample_name'],
        config = config['config']['relocate']['files'][0]

    output:
        config['outdir']['relocate']+"unfiltered/combined.gff"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE/relocate_run.py"

rule relocaTE_post:
    input:
        relocate_gff = config['outdir']['relocate']+"unfiltered/combined.gff",
        te_gff = config['mcc']['relocaTE_ref_TEs'],
        reference_fasta = config['mcc']['reference']

    threads: 1

    conda: config['envs']['processing']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['outdir']['relocate'],
        log = config['args']['log_dir']+"relocaTE.log",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        config = config['config']['relocate']['files'][1]

    output:
        config['out']['relocate']
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE/relocate_post.py"