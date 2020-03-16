
rule setup_reads:
    output:
        config['mcc']['fq1'],
        config['mcc']['fq2']
        
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
        config['mcc']['consensus']
    
    run:
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+config['in']['reference']+" 80 > "+config['mcc']['reference'])
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+config['in']['consensus']+" 80 > "+config['mcc']['consensus'])


rule coverage:
    input:
        config['mcc']['fq1'],
        config['mcc']['fq2'],
        config['mcc']['reference'],
        config['mcc']['consensus']
    
    output:
        config['args']['out']+"/coverage/coverage.log"
    
    run:
        shell("touch "+config['args']['out']+"/coverage/coverage.log")