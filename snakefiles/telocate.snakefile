rule telocate_taxonomy:
    input:
        script = config['args']['mcc_path']+"/install/tools/te-locate/TE_hierarchy.pl",
        ref_gff = config['mcc']['locations'],
        taxonomy = config['mcc']['taxonomy']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log",
        tmp_dir=config['args']['out']+"/tmp"

    conda: config['envs']['processing']

    output:
        config['mcc']['telocate_te_gff']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/telocate_taxonomy.py"

rule telocate_sam:
    input:
        config['mcc']['sam']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['processing']
    
    output:
        config['mcc']['telocate_sam']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/telocate_sam.py"

rule telocate_ref:
    input:
        config['mcc']['reference']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['processing']  

    output:
        config['mcc']['telocate_ref_fasta']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/telocate_ref.py"

rule telocate_run:
    input:
        te_gff = config['mcc']['telocate_te_gff'],
        sam = config['mcc']['telocate_sam'],
        ref =config['mcc']['telocate_ref_fasta'],
        median_insert_size = config['mcc']['median_insert_size']
    
    params:
        run_script = config['args']['mcc_path']+"/install/tools/te-locate/TE_locate.pl",
        out_dir = config['args']['out']+"/results/te-locate/unfiltered/",
        log = config['args']['log_dir']+"te-locate.log"
    
    threads: 1

    conda: config['envs']['te-locate']

    output:
        config['args']['out']+"results/te-locate/unfiltered/te-locate-raw.info"
    
    script:
        config['args']['mcc_path']+"/scripts/telocate/telocate_run.py"

rule telocate_post:
    input:
        telocate_raw = config['args']['out']+"results/te-locate/unfiltered/te-locate-raw.info",
        te_gff = config['mcc']['telocate_te_gff'],
        reference_fasta = config['mcc']['reference']
    
    params:
        out_dir = config['args']['out']+"/results/te-locate/",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']
    
    threads: 1

    conda: config['envs']['processing']

    output:
        config['out']['te-locate']
    
    script:
        config['args']['mcc_path']+"/scripts/telocate/telocate_post.py"