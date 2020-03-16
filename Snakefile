
rule setup_reads:
    output:
        config['args']['out']+"/inputs/"+config['args']['sample_name']+"_1.fastq",
        config['args']['out']+"/inputs/"+config['args']['sample_name']+"_2.fastq"
        
    run:
        if config['args']['fq2'] == "None":
            shell("touch "+config['args']['out']+"/inputs/"+config['args']['sample_name']+"_2.fastq")
        else:
            shell("cp "+config['args']['fq2']+" "+config['args']['out']+"/inputs/"+config['args']['sample_name']+"_2.fastq")
        
        shell("cp "+config['args']['fq1']+" "+config['args']['out']+"/inputs/"+config['args']['sample_name']+"_1.fastq")
        shell("cp "+config['args']['reference']+" "+config['args']['out']+"/inputs/"+config['args']['ref_name']+".fasta")
        shell("cp "+config['args']['consensus']+" "+config['args']['out']+"/inputs/consensusTEs.fasta")


rule fix_line_lengths:
    output:
        config['args']['out']+"/inputs/"+config['args']['ref_name']+".fasta",
        config['args']['out']+"/inputs/consensusTEs.fasta"
    
    run:
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+config['args']['consensus']+" 80 > "+config['args']['out']+"/inputs/consensusTEs.fasta")
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+config['args']['reference']+" 80 > "+config['args']['out']+"/inputs/"+config['args']['ref_name']+".fasta")


rule coverage:
    input:
        config['args']['out']+"/inputs/"+config['args']['sample_name']+"_1.fastq",
        config['args']['out']+"/inputs/"+config['args']['sample_name']+"_2.fastq",
        config['args']['out']+"/inputs/"+config['args']['ref_name']+".fasta",
        config['args']['out']+"/inputs/consensusTEs.fasta"
    
    output:
        config['args']['out']+"/coverage/coverage.log"
    
    run:
        shell("touch "+config['args']['out']+"/coverage/coverage.log")