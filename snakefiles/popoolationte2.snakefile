
rule popoolationTE2_preprocessing:
    input:
        ref_fasta = config['mcc']['popoolationTE_ref_fasta'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['popoolationte2']

    params:
        out_dir = config['outdir']['popoolationte2']+"unfiltered/",
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"popoolationTE2.log",
        jar = config['args']['mcc_path']+"/install/tools/popoolationte2/popte2-v1.10.03.jar"

    output:
        bam = config['outdir']['popoolationte2']+"unfiltered/sorted.bam"

    script:
        config['args']['mcc_path']+"/scripts/popoolationte2/popoolationte2_pre.py"

rule popoolationTE2_run:
    input:
        ref_fasta = config['mcc']['popoolationTE_ref_fasta'],
        taxonomy = config['mcc']['popoolationTE_taxonomy'],
        te_gff = config['mcc']['popoolationTE_gff'],
        bam = config['outdir']['popoolationte2']+"unfiltered/sorted.bam"
    
    threads: 1

    conda: config['envs']['popoolationte2']

    params:
        out_dir = config['outdir']['popoolationte2']+"unfiltered/",
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"popoolationTE2.log",
        jar = config['args']['mcc_path']+"/install/tools/popoolationte2/popte2-v1.10.03.jar",
        config = config['config']['popoolationte2']['files'][0]

    output:
        config['outdir']['popoolationte2']+"unfiltered/teinsertions.txt"

    script:
        config['args']['mcc_path']+"/scripts/popoolationte2/popoolationte2_run.py"


rule popoolationTE2_post:
    input:
        popoolationte2_out = config['outdir']['popoolationte2']+"unfiltered/teinsertions.txt",
        te_gff = config['mcc']['popoolationTE_gff'],
        taxonomy = config['mcc']['popoolationTE_taxonomy'],
        reference_fasta = config['mcc']['popoolationTE_ref_fasta']
    
    threads: 1

    conda: config['envs']['processing']

    params:
        out_dir = config['outdir']['popoolationte2'],
        sample_name=config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        log = config['args']['log_dir']+"popoolationTE2.log",
        config = config['config']['popoolationte2']['files'][1]

    output:
        config['out']['popoolationte2']

    script:
        config['args']['mcc_path']+"/scripts/popoolationte2/popoolationte2_post.py"
