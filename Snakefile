
rule setup_reads:
    output:
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: 1
        
    run:
        # fastq1
        if ".gz" in config['in']['fq1']:
            shell("zcat "+config['in']['fq1']+" > "+config['mcc']['fq1'])
        else:
            shell("cp "+config['in']['fq1']+" "+config['mcc']['fq1'])

        # fastq2
        if config['in']['fq2'] == "None":
            shell("touch "+config['mcc']['fq2'])
        else:
            if ".gz" in config['in']['fq2']:
                shell("zcat "+config['in']['fq2']+" > "+config['mcc']['fq2'])
            else:
                shell("cp "+config['in']['fq2']+" "+config['mcc']['fq2'])



rule fix_line_lengths:
    output:
        config['mcc']['reference'],
        temp(config['mcc']['consensus']),
        temp(config['mcc']['coverage_fasta'])

    threads: 1
    
    run:
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+config['in']['reference']+" 80 > "+config['mcc']['reference'])
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+config['in']['consensus']+" 80 > "+config['mcc']['consensus'])
        if config['in']['coverage_fasta'] == "None":
            shell("touch "+config['mcc']['coverage_fasta'])
        else:
            shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+config['in']['coverage_fasta']+" 80 > "+config['mcc']['coverage_fasta'])


rule make_run_copies:
    output:
        temp(config['mcc']['locations']),
        temp(config['mcc']['taxonomy'])
    run:
        if config['in']['locations'] == "None":
            shell("touch "+config['mcc']['locations'])
        else:
            shell("cp "+config['in']['locations']+" "+config['mcc']['locations'])
        
        if config['in']['taxonomy'] == "None":
            shell("touch "+config['mcc']['taxonomy'])
        else:
            shell("cp "+config['in']['taxonomy']+" "+config['mcc']['taxonomy'])

rule make_reference_te_gff:
    input:
        config['mcc']['reference'],
        config['mcc']['consensus'],
        config['mcc']['locations'],
        config['mcc']['taxonomy']
    
    threads: workflow.cores

    output:
        config['mcc']['masked_fasta'],
        config['mcc']['formatted_ref_tes'],
        config['mcc']['formatted_taxonomy'],
        config['mcc']['formatted_consensus'],
        config['mcc']['ref_te_fasta'],
        config['mcc']['augmented_reference']
    
    script:
        config['args']['mcc_path']+"/scripts/make_ref_te_files.py"

rule coverage:
    input:
        config['mcc']['fq1'],
        config['mcc']['fq2'],
        config['mcc']['reference'],
        config['mcc']['consensus'],
        config['mcc']['coverage_fasta']
    
    output:
        config['args']['out']+"/coverage/coverage.log"

    run:
        shell("touch "+config['args']['out']+"/coverage/coverage.log")    
    # script:
    #     "modules/coverage.py"

rule telocate:
    input:
        config['mcc']['formatted_ref_tes']
    
    output:
        config['args']['out']+"/te-locate/te-locate.log"
    
    run:
        shell("touch "+config['args']['out']+"/te-locate/te-locate.log")