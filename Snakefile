
rule setup_reads:
    input:
        config['in']['fq1'],
        config['in']['fq2']

    params:
        log=config['args']['out']+"/logs/processing.log"

    output:
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: workflow.cores

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"
        
    script:
        config['args']['mcc_path']+"/scripts/setup_reads.py"


rule fix_line_lengths:
    input:
        ref = config['in']['reference'],
        consensus = config['in']['consensus']

    params:
        coverage_fasta = config['in']['coverage_fasta'],
        log=config['args']['out']+"/logs/processing.log"

    output:
        temp(config['mcc']['reference']),
        config['mcc']['consensus'],
        temp(config['mcc']['coverage_fasta'])

    threads: 1

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    script:
        config['args']['mcc_path']+"/scripts/fix_line_lengths.py"


rule make_run_copies:
    params:
        log=config['args']['out']+"/logs/processing.log"

    output:
        temp(config['mcc']['locations']),
        temp(config['mcc']['taxonomy'])

    threads: 1

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    script:
        config['args']['mcc_path']+"/scripts/make_run_copies.py"

rule make_reference_te_files:
    input:
        config['mcc']['reference'],
        config['mcc']['consensus'],
        config['mcc']['locations'],
        config['mcc']['taxonomy']
    
    threads: workflow.cores

    params:
        log=config['args']['out']+"/logs/processing.log"

    output:
        config['mcc']['masked_fasta'],
        config['mcc']['formatted_ref_tes'],
        config['mcc']['formatted_taxonomy'],
        config['mcc']['formatted_consensus'],
        config['mcc']['ref_te_fasta'],
        config['mcc']['augmented_reference'],
        config['mcc']['popoolationTE_taxonomy'],
        config['mcc']['popoolationTE_gff']
    
    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"
    
    script:
        config['args']['mcc_path']+"/scripts/make_reference_te_files.py"

rule index_reference_genome:
    input:
        ref = config['mcc']['augmented_reference']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        config['mcc']['augmented_reference']+".fai",
        config['mcc']['augmented_reference']+".bwt"
    
    script:
        config['args']['mcc_path']+"/scripts/index_reference_genome.py"


rule map_reads:
    input:
        ref = config['mcc']['augmented_reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        idx = config['mcc']['augmented_reference']+".bwt"
    
    params:
        sample=config['args']['sample_name'],
        log=config['args']['out']+"/logs/processing.log"
    
    log: config['args']['out']+"/logs/bwa.log"
    
    threads: workflow.cores

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output: config['mcc']['sam']

    script:
        config['args']['mcc_path']+"/scripts/map_reads.py"

rule sam_to_bam:
    input:
        sam = config['mcc']['sam'],
        ref_idx = config['mcc']['augmented_reference']+".fai"

    params:
        log=config['args']['out']+"/logs/processing.log"
    
    threads: workflow.cores

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        bam = config['mcc']['bam'],
        flagstat = config['summary']['flagstat'],
        tmp_bam = temp(config['mcc']['mcc_files']+config['args']['run_id']+".tmp.bam")
    
    script:
        config['args']['mcc_path']+"/scripts/sam_to_bam.py"
        

rule make_ref_te_bed:
    input:
        config['mcc']['formatted_ref_tes']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        config['mcc']['ref_tes_bed']

    script:
        config['args']['mcc_path']+"/scripts/make_ref_te_bed.py"

rule telocate_taxonomy:
    input:
        script = config['args']['mcc_path']+"/install/tools/te-locate/TE_hierarchy.pl",
        ref_gff = config['mcc']['formatted_ref_tes'],
        taxonomy = config['mcc']['formatted_taxonomy']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        config['mcc']['telocate_te_gff']
    
    script:
        config['args']['mcc_path']+"/scripts/telocate_taxonomy.py"


rule median_insert_size:
    input:
        config['mcc']['sam']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        config['summary']['median_insert_size']
    
    script:
        config['args']['mcc_path']+"/scripts/median_insert_size.py"

rule telocate_sam:
    input:
        config['mcc']['sam']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"
    
    output:
        config['mcc']['telocate_sam']

    script:
        config['args']['mcc_path']+"/scripts/telocate_sam.py"

rule telocate_ref:
    input:
        config['mcc']['augmented_reference']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"    

    output:
        config['mcc']['telocate_ref_fasta']
    
    script:
        config['args']['mcc_path']+"/scripts/telocate_ref.py"

rule reference_2bit:
    input:
        config['mcc']['augmented_reference']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"
    
    output:
        config['mcc']['ref_2bit']
    
    script:
        config['args']['mcc_path']+"/scripts/reference_2bit.py"

rule relocaTE_consensus:
    input:
        config['mcc']['formatted_consensus']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        config['mcc']['relocaTE_consensus']

    script:
        config['args']['mcc_path']+"/scripts/relocaTE_consensus.py"
    
rule relocaTE_ref_gff:
    input:
        config['mcc']['formatted_ref_tes'],
        config['mcc']['formatted_taxonomy']

    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        config['mcc']['relocaTE_ref_TEs']

    script:
        config['args']['mcc_path']+"/scripts/relocaTE_ref_gff.py" 


rule popoolationTE_ref_fasta:
    input:
        config['mcc']['masked_fasta'],
        config['mcc']['formatted_consensus'],
        config['mcc']['ref_te_fasta']
    
    threads: 1

    params:
        log=config['args']['out']+"/logs/processing.log"

    conda: config['args']['mcc_path']+"/envs/mcc_processing.yml"

    output:
        config['mcc']['popoolationTE_ref_fasta']
    
    script:
        config['args']['mcc_path']+"/scripts/popoolationTE_ref_fasta.py" 


rule coverage:
    input:
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        ref = config['mcc']['reference'],
        consensus = config['mcc']['consensus'],
        coverage_fa = config['mcc']['coverage_fasta']
    
    params: 
        sample=config['args']['sample_name'],
        log=config['args']['out']+"/logs/coverage.log"

    threads: workflow.cores

    conda: config['args']['mcc_path']+"/envs/mcc_coverage.yml"

    output:
        config['args']['out']+"/results/coverage/te_depth.csv"

    script:
        config['args']['mcc_path']+"/modules/coverage/coverage.py"   


rule telocate:
    input:
        config['mcc']['telocate_te_gff'],
        config['summary']['median_insert_size'],
        config['mcc']['telocate_sam'],
        config['mcc']['telocate_ref_fasta'],
    
    threads: 1

    output:
        config['args']['out']+"/results/te-locate/te-locate.log"
    
    run:
        shell("touch "+config['args']['out']+"/results/te-locate/te-locate.log")


rule retroseq:
    input:
        config['mcc']['formatted_consensus'],
        config['mcc']['bam'],
        config['mcc']['augmented_reference'],
        config['mcc']['ref_tes_bed'],
        config['mcc']['formatted_taxonomy']

    threads: 1
    
    output:
        config['args']['out']+"/results/retroseq/retroseq.log"
    
    run:
        shell("touch "+config['args']['out']+"/results/retroseq/retroseq.log")



rule run_temp:
    input:
        config['args']['mcc_path']+"/config/TEMP/temp_run.py",
        bam = config['mcc']['bam'],
        twobit = config['mcc']['ref_2bit'],
        consensus = config['mcc']['formatted_consensus'],
        ref_te_bed = config['mcc']['ref_tes_bed'],
        taxonomy = config['mcc']['formatted_taxonomy'],
        median_insert_size = config['summary']['median_insert_size']
        
    
    conda: config['args']['mcc_path']+"/envs/mcc_temp.yml"

    params:
        log = config['args']['out']+"/logs/TEMP.log",
        scripts_dir = config['args']['mcc_path']+"/install/tools/temp/scripts/",
        out_dir = config['args']['out']+"/results/TEMP/unfiltered/",
        sample_name = config['args']['sample_name']

    threads: workflow.cores

    output:
        config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary"
    
    script:
        config['args']['mcc_path']+"/modules/TEMP/temp_run.py"

rule process_temp:
    input:
        config['args']['mcc_path']+"/config/TEMP/temp_post.py",
        insert_summary = config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        absence_summary = config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary",
        te_gff = config['mcc']['telocate_te_gff']
    
    conda: config['args']['mcc_path']+"/envs/mcc_temp.yml"

    params:
        log = config['args']['out']+"/logs/TEMP.log",
        out_dir = config['args']['out']+"/results/TEMP/",
        sample_name = config['args']['sample_name']

    threads: workflow.cores

    output:
        config['args']['out']+"/results/TEMP/"+config['args']['sample_name']+"_temp_raw.bed",
        config['args']['out']+"/results/TEMP/"+config['args']['sample_name']+"_temp_redundant.bed",
        config['args']['out']+"/results/TEMP/"+config['args']['sample_name']+"_temp_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/modules/TEMP/temp_post.py"


rule relocaTE:
    input:
        config['mcc']['relocaTE_consensus'],
        config['mcc']['relocaTE_ref_TEs'],
        config['mcc']['augmented_reference'],
        config['mcc']['fq1'],
        config['mcc']['fq2']

    threads: workflow.cores

    output:
        config['args']['out']+"/results/relocaTE/relocaTE.log"
    
    run:
        shell("touch "+output[0])


rule ngs_te_mapper:
    input:
        consensus_fasta = config['mcc']['formatted_consensus'],
        reference_fasta = config['mcc']['augmented_reference'],
        fastq1 = config['mcc']['fq1'],
        fastq2 = config['mcc']['fq2']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/ngs_te_mapper/",
        script_dir = config['args']['mcc_path']+"/install/tools/ngs_te_mapper/sourceCode/",
        log = config['args']['out']+"/logs/ngs_te_mapper.log",
        sample_name = config['args']['sample_name']
    
    threads: workflow.cores

    conda: config['args']['mcc_path']+"/envs/mcc_ngs_te_mapper.yml"

    output:
        config['args']['out']+"/results/ngs_te_mapper/"+config['args']['sample_name']+"_ngs_te_mapper_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/modules/ngs_te_mapper/ngs_te_mapper.py"
        

rule popoolationTE:
    input:
        config['mcc']['popoolationTE_ref_fasta'],
        config['mcc']['popoolationTE_taxonomy'],
        config['mcc']['popoolationTE_gff'],
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: workflow.cores

    output:
        config['args']['out']+"/results/popoolationTE/popoolationTE.log"

    run:
        shell("touch "+output[0])

