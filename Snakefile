include: config['args']['mcc_path']+"/snakefiles/temp.snakefile"
include: config['args']['mcc_path']+"/snakefiles/relocate.snakefile"
include: config['args']['mcc_path']+"/snakefiles/relocate2.snakefile"
include: config['args']['mcc_path']+"/snakefiles/popoolationte.snakefile"
include: config['args']['mcc_path']+"/snakefiles/popoolationte2.snakefile"
include: config['args']['mcc_path']+"/snakefiles/telocate.snakefile"
include: config['args']['mcc_path']+"/snakefiles/ngs_te_mapper.snakefile"
include: config['args']['mcc_path']+"/snakefiles/retroseq.snakefile"
include: config['args']['mcc_path']+"/snakefiles/tepid.snakefile"

rule setup_reads:
    input:
        fq1 = config['in']['fq1']
        
    params:
        fq2 = config['in']['fq2'],
        methods = config['args']['methods'],
        out = config['args']['out'],
        run_id = config['args']['run_id'],
        log = config['args']['log_dir']+"trimgalore.log"

    output:
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']
        
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/setup_reads.py"


rule make_coverage_fasta:
    params:
        coverage_fasta = config['in']['coverage_fasta']
    
    threads: 1

    conda: config['envs']['mcc_processing']

    output:
        coverage_fasta = config['mcc']['coverage_fasta']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_coverage_fasta.py"


rule make_reference_fasta:
    input:
        ref = config['in']['reference']
    
    threads: 1

    params: 
        log = config['args']['log_dir']+"processing.log",
        augment = config['args']['augment_fasta'],
        mcc_out = config['args']['out'],
        run_id = config['args']['run_id']
    
    output: 
        ref = config['mcc']['unaugmented_reference'],
        aug_ref = config['mcc']['reference']

    conda: 
        config['envs']['mcc_processing']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_reference_fasta.py"

rule make_consensus_fasta:
    input:
        consensus = config['in']['consensus']
    
    params:
        mcc_out = config['args']['out'],
        run_id = config['args']['run_id']

    threads: 1

    conda:
        config['envs']['mcc_processing']
    
    output:
        consensus = config['mcc']['consensus']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_consensus_fasta.py"

rule make_te_annotations:
    input:
        ref = config['mcc']['unaugmented_reference'],
        consensus = config['mcc']['consensus']

    params:
        te_gff = config['in']['locations'],
        taxonomy = config['in']['taxonomy'],
        mcc_out = config['args']['out'],
        run_id = config['args']['run_id'],
        log = config['args']['log_dir']+"processing.log",
        chromosomes = config['args']['chromosomes'],
        augment = config['args']['augment_fasta']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']

    output:
        te_gff = config['mcc']['unaugmented_locations'],
        aug_te_gff = config['mcc']['locations'],
        taxonomy = config['mcc']['unaugmented_taxonomy'],
        aug_taxonomy = config['mcc']['taxonomy']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_te_annotations.py"


rule mask_reference_fasta:
    input:
        ref_fasta = config['mcc']['unaugmented_reference'],
        te_gff = config['mcc']['unaugmented_locations']
    
    params:
        mcc_out = config['args']['out'],
        run_id = config['args']['run_id'],
        log = config['args']['log_dir']+"processing.log",
        chromosomes = config['args']['chromosomes'],
        augment = config['args']['augment_fasta']
    
    threads: 1

    conda: config['envs']['mcc_processing']

    output:
        masked_ref = config['mcc']['masked_fasta']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/mask_reference_fasta.py"

rule make_ref_te_fasta:
    input:
        ref_fasta = config['mcc']['unaugmented_reference'],
        te_gff = config['mcc']['unaugmented_locations']
    
    params:
        mcc_out = config['args']['out'],
        run_id = config['args']['run_id'],
        log = config['args']['log_dir']+"processing.log",
        chromosomes = config['args']['chromosomes'],
        augment = config['args']['augment_fasta']
    
    threads: 1

    conda: config['envs']['mcc_processing']

    output:
        ref_te_fasta = config['mcc']['ref_te_fasta']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_ref_te_fasta.py"


rule index_reference_genome:
    input:
        ref = config['mcc']['reference']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['reference']+".fai",
        config['mcc']['reference']+".bwt"
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/index_reference_genome.py"


rule map_reads:
    input:
        ref = config['mcc']['reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        idx = config['mcc']['reference']+".bwt"
    
    params:
        sample=config['args']['sample_name'],
        log=config['args']['log_dir']+"processing.log"
    
    log: config['args']['log_dir']+"bwa.log"
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']

    output: config['mcc']['sam']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/map_reads.py"

rule sam_to_bam:
    input:
        sam = config['mcc']['sam'],
        ref_idx = config['mcc']['reference']+".fai"

    params:
        log=config['args']['log_dir']+"processing.log"
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']

    output:
        bam = config['mcc']['bam'],
        flagstat = config['mcc']['flagstat'],
        tmp_bam = temp(config['mcc']['mcc_files']+config['args']['run_id']+".tmp.bam")
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/sam_to_bam.py"
        

rule make_ref_te_bed:
    input:
        config['mcc']['locations']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['ref_tes_bed']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_ref_te_bed.py"

rule median_insert_size:
    input:
        config['mcc']['sam']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log",
        fq2 = config['in']['fq2']

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['median_insert_size']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/median_insert_size.py"

rule reference_2bit:
    input:
        config['mcc']['reference']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']
    
    output:
        config['mcc']['ref_2bit']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/reference_2bit.py"

rule repeatmask:
    input:
        reference = config['mcc']['reference'],
        te_seqs = config['mcc']['consensus']
    
    params:
        ref_name = config['args']['ref_name'],
        out_dir = config['args']['out'],
        log=config['args']['log_dir']+"processing.log"
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']

    output:
        rm_out = config['mcc']['repeatmasker_out']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/repeatmask.py"


rule coverage:
    input:
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        ref = config['mcc']['unaugmented_reference'],
        consensus = config['mcc']['consensus'],
        coverage_fa = config['mcc']['coverage_fasta']
    
    params: 
        sample=config['args']['sample_name'],
        log=config['args']['log_dir']+"coverage.log"

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['coverage']

    output:
        config['args']['out']+"/results/coverage/te_depth.csv"

    script:
        config['args']['mcc_path']+"/scripts/coverage/coverage.py"   

rule summary_report:
    input:
        out_files = config['args']['out_files'].split(","),
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']

    params:
        commit = config['args']['commit'],
        ref = config['mcc']['reference'],
        consensus = config['mcc']['consensus'],
        taxonomy = config['mcc']['taxonomy'],
        bam = config['mcc']['bam'],
        flagstat = config['mcc']['flagstat'],
        median_insert_size = config['mcc']['median_insert_size'],
        methods = config['args']['methods'].split(","),
        results_dir = config['args']['out']+"/results/",
        sample_name = config['args']['sample_name'],
        command = config['args']['full_command'],
        execution_dir = config['args']['call_directory'],
        time = config['args']['time'],
        raw_fq1 = config['in']['fq1'],
        raw_fq2 = config['in']['fq2'],
        chromosomes = config['args']['chromosomes'],
        out_dir = config['args']['out']+"/results/summary/"


    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']

    output:
        summary_report = config['args']['out']+"/results/summary/data/run/summary_report.txt",
        html_summary_report = config['args']['out']+"/results/summary/summary.html"
    
    script:
        config['args']['mcc_path']+"/scripts/summary/summary_report.py"
        