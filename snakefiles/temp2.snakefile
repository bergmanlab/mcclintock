rule run_temp2:
    input:
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        bam = config['mcc']['bam'],
        reference = config['mcc']['reference'],
        twobit = config['mcc']['ref_2bit'],
        consensus = config['mcc']['consensus'],
        ref_te_bed = config['mcc']['ref_tes_bed'],
        taxonomy = config['mcc']['taxonomy'],
        median_insert_size = config['mcc']['median_insert_size']
        
    
    conda: config['envs']['temp2']

    params:
        log = config['args']['log_dir']+"temp2.log",
        script_dir = config['args']['mcc_path']+"/install/tools/temp2/",
        out_dir = config['outdir']['temp2']+"unfiltered/",
        sample_name = config['args']['sample_name'],
        config = config['config']['temp2']['files'][0],
        status_log = config['status']['temp2']

    threads: config['args']['max_threads_per_rule']

    output:
        config['outdir']['temp2']+"unfiltered/"+config['args']['sample_name']+".insertion.bed",
        config['outdir']['temp2']+"unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary"
    
    script:
        config['args']['mcc_path']+"/scripts/temp2/temp2_run.py"

rule process_temp2:
    input:
        insert_bed = config['outdir']['temp2']+"unfiltered/"+config['args']['sample_name']+".insertion.bed",
        absence_summary = config['outdir']['temp2']+"unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary",
        te_gff = config['mcc']['telocate_te_gff'],
        reference_fasta = config['mcc']['reference']
    
    conda: config['envs']['processing']

    params:
        log = config['args']['log_dir']+"temp2.log",
        out_dir = config['outdir']['temp2'],
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        config = config['config']['temp2']['files'][1],
        status_log = config['status']['temp2']

    threads: 1

    output:
        config['out']['temp2']
    
    script:
        config['args']['mcc_path']+"/scripts/temp2/temp2_post.py"