rule run_temp:
    input:
        bam = config['mcc']['bam'],
        twobit = config['mcc']['ref_2bit'],
        consensus = config['mcc']['consensus'],
        ref_te_bed = config['mcc']['ref_tes_bed'],
        taxonomy = config['mcc']['taxonomy'],
        median_insert_size = config['mcc']['median_insert_size']
        
    
    conda: config['envs']['temp']

    params:
        log = config['args']['log_dir']+"TEMP.log",
        scripts_dir = config['args']['mcc_path']+"/install/tools/temp/scripts/",
        out_dir = config['outdir']['temp']+"unfiltered/",
        sample_name = config['args']['sample_name'],
        config = config['config']['temp']['files'][0],
        status_log = config['status']['temp']

    threads: config['args']['max_threads_per_rule']

    output:
        config['outdir']['temp']+"unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        config['outdir']['temp']+"unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary"
    
    script:
        config['args']['mcc_path']+"/scripts/TEMP/temp_run.py"

rule process_temp:
    input:
        insert_summary = config['outdir']['temp']+"unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        absence_summary = config['outdir']['temp']+"unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary",
        te_gff = config['mcc']['telocate_te_gff'],
        reference_fasta = config['mcc']['reference']
    
    conda: config['envs']['processing']

    params:
        log = config['args']['log_dir']+"TEMP.log",
        out_dir = config['outdir']['temp'],
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes'],
        config = config['config']['temp']['files'][1],
        status_log = config['status']['temp']

    threads: 1

    output:
        config['out']['temp']
    
    script:
        config['args']['mcc_path']+"/scripts/TEMP/temp_post.py"