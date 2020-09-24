rule run_temp:
    input:
        config['args']['mcc_path']+"/config/TEMP/temp_run.py",
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
        out_dir = config['args']['out']+"/results/TEMP/unfiltered/",
        sample_name = config['args']['sample_name']

    threads: config['args']['max_threads_per_rule']

    output:
        config['args']['out']+"results/TEMP/unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        config['args']['out']+"results/TEMP/unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary"
    
    script:
        config['args']['mcc_path']+"/scripts/TEMP/temp_run.py"

rule process_temp:
    input:
        config['args']['mcc_path']+"/config/TEMP/temp_post.py",
        insert_summary = config['args']['out']+"results/TEMP/unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        absence_summary = config['args']['out']+"results/TEMP/unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary",
        te_gff = config['mcc']['telocate_te_gff']
    
    conda: config['envs']['temp']

    params:
        log = config['args']['log_dir']+"TEMP.log",
        out_dir = config['args']['out']+"/results/TEMP/",
        sample_name = config['args']['sample_name'],
        chromosomes = config['args']['chromosomes']

    threads: 1

    output:
        config['out']['temp']
    
    script:
        config['args']['mcc_path']+"/scripts/TEMP/temp_post.py"