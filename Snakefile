
rule setup_reads:
    input:
        config['in']['fq1'],
        

    params:
        fq2 = config['in']['fq2'],
        log = config['args']['log_dir']+"trimgalore.log"

    output:
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']
        
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/setup_reads.py"


rule fix_line_lengths:
    input:
        ref = config['in']['reference'],
        consensus = config['in']['consensus']

    params:
        coverage_fasta = config['in']['coverage_fasta'],
        log=config['args']['log_dir']+"processing.log"

    output:
        config['mcc']['reference'],
        config['mcc']['consensus'],
        config['mcc']['coverage_fasta']

    threads: 1

    conda: config['envs']['mcc_processing']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/fix_line_lengths.py"


rule make_run_copies:
    params:
        log=config['args']['log_dir']+"processing.log"

    output:
        config['mcc']['locations'],
        config['mcc']['taxonomy']

    threads: 1

    conda: config['envs']['mcc_processing']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_run_copies.py"

rule make_reference_te_files:
    input:
        config['mcc']['reference'],
        config['mcc']['consensus'],
        config['mcc']['locations'],
        config['mcc']['taxonomy']
    
    threads: config['args']['max_threads_per_rule']

    params:
        log=config['args']['log_dir']+"processing.log",
        raw_ref_tes = config['in']['locations'],
        out_dir = config['args']['out'],
        run_id = config['args']['run_id']

    output:
        config['mcc']['masked_fasta'],
        config['mcc']['formatted_ref_tes'],
        config['mcc']['formatted_taxonomy'],
        config['mcc']['formatted_consensus'],
        config['mcc']['ref_te_fasta'],
        config['mcc']['augmented_reference'],
        config['mcc']['popoolationTE_taxonomy'],
        config['mcc']['popoolationTE_gff']
    
    conda: config['envs']['mcc_processing']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_reference_te_files.py"

rule index_reference_genome:
    input:
        ref = config['mcc']['augmented_reference']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['augmented_reference']+".fai",
        config['mcc']['augmented_reference']+".bwt"
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/index_reference_genome.py"


rule map_reads:
    input:
        ref = config['mcc']['augmented_reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        idx = config['mcc']['augmented_reference']+".bwt"
    
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
        ref_idx = config['mcc']['augmented_reference']+".fai"

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
        config['mcc']['formatted_ref_tes']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['ref_tes_bed']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/make_ref_te_bed.py"

rule telocate_taxonomy:
    input:
        script = config['args']['mcc_path']+"/install/tools/te-locate/TE_hierarchy.pl",
        ref_gff = config['mcc']['formatted_ref_tes'],
        taxonomy = config['mcc']['formatted_taxonomy']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['telocate_te_gff']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/telocate_taxonomy.py"


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

rule telocate_sam:
    input:
        config['mcc']['sam']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']
    
    output:
        config['mcc']['telocate_sam']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/telocate_sam.py"

rule telocate_ref:
    input:
        config['mcc']['augmented_reference']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']  

    output:
        config['mcc']['telocate_ref_fasta']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/telocate_ref.py"

rule reference_2bit:
    input:
        config['mcc']['augmented_reference']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']
    
    output:
        config['mcc']['ref_2bit']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/reference_2bit.py"

rule relocaTE_consensus:
    input:
        config['mcc']['formatted_consensus']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['relocaTE_consensus']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/relocaTE_consensus.py"
    
rule relocaTE_ref_gff:
    input:
        config['mcc']['formatted_ref_tes'],
        config['mcc']['formatted_taxonomy']

    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['relocaTE_ref_TEs']

    script:
        config['args']['mcc_path']+"/scripts/preprocessing/relocaTE_ref_gff.py" 


rule popoolationTE_ref_fasta:
    input:
        config['mcc']['masked_fasta'],
        config['mcc']['formatted_consensus'],
        config['mcc']['ref_te_fasta']
    
    threads: 1

    params:
        log=config['args']['log_dir']+"processing.log"

    conda: config['envs']['mcc_processing']

    output:
        config['mcc']['popoolationTE_ref_fasta']
    
    script:
        config['args']['mcc_path']+"/scripts/preprocessing/popoolationTE_ref_fasta.py"


rule repeatmask:
    input:
        reference = config['mcc']['augmented_reference'],
        te_seqs = config['mcc']['formatted_consensus']
    
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
        ref = config['mcc']['reference'],
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


rule run_temp:
    input:
        config['args']['mcc_path']+"/config/TEMP/temp_run.py",
        bam = config['mcc']['bam'],
        twobit = config['mcc']['ref_2bit'],
        consensus = config['mcc']['formatted_consensus'],
        ref_te_bed = config['mcc']['ref_tes_bed'],
        taxonomy = config['mcc']['formatted_taxonomy'],
        median_insert_size = config['mcc']['median_insert_size']
        
    
    conda: config['envs']['temp']

    params:
        log = config['args']['log_dir']+"TEMP.log",
        scripts_dir = config['args']['mcc_path']+"/install/tools/temp/scripts/",
        out_dir = config['args']['out']+"/results/TEMP/unfiltered/",
        sample_name = config['args']['sample_name']

    threads: config['args']['max_threads_per_rule']

    output:
        config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary"
    
    script:
        config['args']['mcc_path']+"/scripts/TEMP/temp_run.py"

rule process_temp:
    input:
        config['args']['mcc_path']+"/config/TEMP/temp_post.py",
        insert_summary = config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".insertion.refined.bp.summary",
        absence_summary = config['args']['out']+"/results/TEMP/unfiltered/"+config['args']['sample_name']+".absence.refined.bp.summary",
        te_gff = config['mcc']['telocate_te_gff']
    
    conda: config['envs']['temp']

    params:
        log = config['args']['log_dir']+"TEMP.log",
        out_dir = config['args']['out']+"/results/TEMP/",
        sample_name = config['args']['sample_name']

    threads: 1

    output:
        config['args']['out']+"/results/TEMP/"+config['args']['sample_name']+"_temp_raw.bed",
        config['args']['out']+"/results/TEMP/"+config['args']['sample_name']+"_temp_redundant.bed",
        config['args']['out']+"/results/TEMP/"+config['args']['sample_name']+"_temp_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/TEMP/temp_post.py"


rule relocaTE_run:
    input:
        config['args']['mcc_path']+"/config/relocate/relocate_run.py",
        consensus_fasta = config['mcc']['relocaTE_consensus'],
        te_gff = config['mcc']['relocaTE_ref_TEs'],
        reference_fasta = config['mcc']['augmented_reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']

    threads: 1

    conda: config['envs']['relocate']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/relocaTE/unfiltered/",
        log = config['args']['log_dir']+"relocaTE.log",
        script_dir = config['args']['mcc_path']+"/install/tools/relocate/scripts/",

    output:
        config['args']['out']+"/results/relocaTE/unfiltered/combined.gff"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE/relocate_run.py"

rule relocaTE_post:
    input:
        config['args']['mcc_path']+"/config/relocate/relocate_post.py",
        relocate_gff = config['args']['out']+"/results/relocaTE/unfiltered/combined.gff",
        te_gff = config['mcc']['relocaTE_ref_TEs'],

    threads: 1

    conda: config['envs']['relocate']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/relocaTE/",
        log = config['args']['log_dir']+"relocaTE.log",
        sample_name = config['args']['sample_name']

    output:
        config['args']['out']+"/results/relocaTE/"+config['args']['sample_name']+"_relocate_redundant.bed",
        config['args']['out']+"/results/relocaTE/"+config['args']['sample_name']+"_relocate_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE/relocate_post.py"


rule relocaTE2_run:
    input:
        config['args']['mcc_path']+"/config/relocate2/relocate2_run.py",
        reference = config['mcc']['augmented_reference'],
        te_seqs = config['mcc']['formatted_consensus'],
        rm_out = config['mcc']['repeatmasker_out'],
        fastq1 = config['mcc']['fq1'],
        fastq2 = config['mcc']['fq2'],
        bam = config['mcc']['bam'],
        median_insert_size = config['mcc']['median_insert_size']

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['relocate2']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/relocaTE2/unfiltered/",
        log = config['args']['log_dir']+"relocaTE2.log"
    
    output:
        config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_nonref_insert.gff",
        config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_ref_insert.gff"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE2/relocate2_run.py"



rule relocaTE2_post:
    input:
        config['args']['mcc_path']+"/config/relocate2/relocate2_post.py",
        rm_out = config['mcc']['repeatmasker_out'],
        nonref_gff = config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_nonref_insert.gff",
        ref_gff = config['args']['out']+"/results/relocaTE2/unfiltered/repeat/results/ALL.all_ref_insert.gff"

    threads: 1

    conda: config['envs']['relocate2']

    params:
        out_dir = config['args']['out']+"/results/relocaTE2/",
        log = config['args']['log_dir']+"relocaTE2.log",
        sample_name = config['args']['sample_name']
    
    output:
        config['args']['out']+"/results/relocaTE2/"+config['args']['sample_name']+"_relocate2_redundant.bed",
        config['args']['out']+"/results/relocaTE2/"+config['args']['sample_name']+"_relocate2_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/relocaTE2/relocate2_post.py"



rule ngs_te_mapper_run:
    input:
        config = config['args']['mcc_path']+"/config/ngs_te_mapper/ngs_te_mapper_run.py",
        consensus_fasta = config['mcc']['formatted_consensus'],
        reference_fasta = config['mcc']['augmented_reference'],
        fastq1 = config['mcc']['fq1'],
        fastq2 = config['mcc']['fq2']

    params:
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/ngs_te_mapper/unfiltered/",
        script_dir = config['args']['mcc_path']+"/install/tools/ngs_te_mapper/sourceCode/",
        log = config['args']['log_dir']+"ngs_te_mapper.log",
        sample_name = config['args']['sample_name']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['ngs_te_mapper']

    output:
        config['args']['out']+"/results/ngs_te_mapper/unfiltered/bed_tsd/"+config['args']['sample_name']+"_insertions.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/ngs_te_mapper/ngs_te_mapper_run.py"

rule ngs_te_mapper_post:
    input:
        config = config['args']['mcc_path']+"/config/ngs_te_mapper/ngs_te_mapper_post.py",
        raw_bed = config['args']['out']+"/results/ngs_te_mapper/unfiltered/bed_tsd/"+config['args']['sample_name']+"_insertions.bed"

    params:
        out_dir = config['args']['out']+"/results/ngs_te_mapper/",
        log = config['args']['log_dir']+"ngs_te_mapper.log",
        sample_name = config['args']['sample_name']
    
    threads: 1

    conda: config['envs']['ngs_te_mapper']

    output:
        config['args']['out']+"/results/ngs_te_mapper/"+config['args']['sample_name']+"_ngs_te_mapper_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/ngs_te_mapper/ngs_te_mapper_post.py"
        


rule telocate_run:
    input:
        te_gff = config['mcc']['telocate_te_gff'],
        sam = config['mcc']['telocate_sam'],
        ref =config['mcc']['telocate_ref_fasta'],
        median_insert_size = config['mcc']['median_insert_size']
    
    params:
        run_script = config['args']['mcc_path']+"/install/tools/te-locate/TE_locate.pl",
        max_mem = config['args']['mem'],
        out_dir = config['args']['out']+"/results/te-locate/unfiltered/",
        log = config['args']['log_dir']+"te-locate.log"
    
    threads: 1

    conda: config['envs']['te-locate']

    output:
        config['args']['out']+"/results/te-locate/unfiltered/te-locate-raw.info"
    
    script:
        config['args']['mcc_path']+"/scripts/telocate/telocate_run.py"

rule telocate_post:
    input:
        telocate_raw = config['args']['out']+"/results/te-locate/unfiltered/te-locate-raw.info",
        te_gff = config['mcc']['telocate_te_gff']
    
    params:
        out_dir = config['args']['out']+"/results/te-locate/",
        sample_name = config['args']['sample_name']
    
    threads: 1

    conda: config['envs']['te-locate']

    output:
        config['args']['out']+"/results/te-locate/"+config['args']['sample_name']+"_telocate_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/telocate/telocate_post.py"


rule retroseq_run:
    input:
        consensus_fasta = config['mcc']['formatted_consensus'],
        bam = config['mcc']['bam'],
        ref_fasta = config['mcc']['augmented_reference'],
        ref_te_bed = config['mcc']['ref_tes_bed'],
        taxonomy = config['mcc']['formatted_taxonomy']

    threads: 1

    conda: config['envs']['retroseq']

    params:
        script_dir = config['args']['mcc_path']+"/install/tools/retroseq/bin/",
        out_dir = config['args']['out']+"/results/retroseq/unfiltered/",
        ref_name=config['args']['ref_name'],
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"retroseq.log"
    
    output:
        config['args']['out']+"/results/retroseq/unfiltered/"+config['args']['sample_name']+".call.PE.vcf"
    
    script:
        config['args']['mcc_path']+"/scripts/retroseq/retroseq_run.py"

rule retroseq_post:
    input:
        retroseq_out = config['args']['out']+"/results/retroseq/unfiltered/"+config['args']['sample_name']+".call.PE.vcf"

    threads: 1

    conda: config['envs']['retroseq']

    params:
        out_dir = config['args']['out']+"/results/retroseq/",
        ref_name=config['args']['ref_name'],
        sample_name=config['args']['sample_name']
    
    output:
        config['args']['out']+"/results/retroseq/"+config['args']['sample_name']+"_retroseq_nonredundant.bed"
    
    script:
        config['args']['mcc_path']+"/scripts/retroseq/retroseq_post.py"


rule popoolationTE_preprocessing:
    input:
        ref_fasta = config['mcc']['popoolationTE_ref_fasta'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2']
    
    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['popoolationte']

    params:
        out_dir = config['args']['out']+"/results/popoolationTE/unfiltered/",
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"popoolationTE.log",
        script_dir = config['args']['mcc_path']+"/install/tools/popoolationte/"

    output:
        config['args']['out']+"/results/popoolationTE/unfiltered/reads1.fastq",
        config['args']['out']+"/results/popoolationTE/unfiltered/reads2.fastq",
        config['args']['out']+"/results/popoolationTE/unfiltered/combined.sorted.sam"

    script:
        config['args']['mcc_path']+"/scripts/popoolationte/popoolationte_pre.py"

rule popoolationTE_run:
    input:
        ref_fasta = config['mcc']['popoolationTE_ref_fasta'],
        taxonomy = config['mcc']['popoolationTE_taxonomy'],
        te_gff = config['mcc']['popoolationTE_gff'],
        fq1 = config['args']['out']+"/results/popoolationTE/unfiltered/reads1.fastq",
        fq2 = config['args']['out']+"/results/popoolationTE/unfiltered/reads2.fastq",
        sam = config['args']['out']+"/results/popoolationTE/unfiltered/combined.sorted.sam"
    
    threads: 1

    conda: config['envs']['popoolationte']

    params:
        out_dir = config['args']['out']+"/results/popoolationTE/unfiltered/",
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"popoolationTE.log",
        script_dir = config['args']['mcc_path']+"/install/tools/popoolationte/"

    output:
        config['args']['out']+"/results/popoolationTE/unfiltered/te-poly-filtered.txt"

    script:
        config['args']['mcc_path']+"/scripts/popoolationte/popoolationte_run.py"

rule popoolationTE_post:
    input:
        popoolationte_out = config['args']['out']+"/results/popoolationTE/unfiltered/te-poly-filtered.txt"
    
    threads: 1

    conda: config['envs']['popoolationte']

    params:
        out_dir = config['args']['out']+"/results/popoolationTE/",
        sample_name=config['args']['sample_name'],
        log = config['args']['log_dir']+"popoolationTE.log"

    output:
        config['args']['out']+"/results/popoolationTE/"+config['args']['sample_name']+"_popoolationte_nonredundant.bed"

    script:
        config['args']['mcc_path']+"/scripts/popoolationte/popoolationte_post.py"


rule summary_report:
    input:
        out_files = config['args']['out_files'].split(","),
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        ref = config['mcc']['reference'],
        bam = config['mcc']['bam'],
        flagstat = config['mcc']['flagstat'],
        median_insert_size = config['mcc']['median_insert_size'],
        taxonomy = config['mcc']['formatted_taxonomy']


    params:
        methods = config['args']['methods'].split(","),
        results_dir = config['args']['out']+"/results/",
        sample_name = config['args']['sample_name'],
        command = config['args']['full_command'],
        execution_dir = config['args']['call_directory'],
        time = config['args']['time'],
        raw_fq2 = config['in']['fq2'],
        out_dir = config['args']['out']+"/results/summary/"


    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['mcc_processing']

    output:
        te_summary = config['args']['out']+"/results/summary/te_summary.txt",
        te_csv = config['args']['out']+"/results/summary/te_summary.csv",
        summary_report = config['args']['out']+"/results/summary/summary_report.txt"
    
    script:
        config['args']['mcc_path']+"/scripts/summary/summary_report.py"
        