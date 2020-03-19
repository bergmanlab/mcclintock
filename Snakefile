
rule setup_reads:
    input:
        config['in']['fq1'],
        config['in']['fq2']
    output:
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: workflow.cores
        
    script:
        config['args']['mcc_path']+"/scripts/trimgalore.py"



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

    threads: 1

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

rule index_reference_genome:
    input:
        config['mcc']['augmented_reference']
    
    threads: 1

    output:
        config['mcc']['augmented_reference']+".fai",
        config['mcc']['augmented_reference']+".bwt"
    
    run:
        shell("samtools faidx "+config['mcc']['augmented_reference'])
        shell("bwa index "+config['mcc']['augmented_reference'])

rule map_reads:
    input:
        ref = config['mcc']['augmented_reference'],
        fq1 = config['mcc']['fq1'],
        fq2 = config['mcc']['fq2'],
        idx = config['mcc']['augmented_reference']+".bwt"
    
    params:
        sample=config['args']['sample_name']
    
    log: config['args']['out']+"/logs/bwa.log"
    
    threads: workflow.cores

    output: config['mcc']['sam']

    run:
        if config['in']['fq2'] == 'None':
            if eval(config['args']['save_comments']):
                shell("bwa mem -C -t "+str(threads)+" -R '@RG\\tID:"+params.sample+"\\tSM:"+params.sample+"' "+input.ref+" "+input.fq1+" > "+output[0]+" 2> "+log[0])
            else:
                shell("bwa mem -t "+str(threads)+" -R '@RG\\tID:"+params.sample+"\\tSM:"+params.sample+"' "+input.ref+" "+input.fq1+" > "+output[0]+" 2> "+log[0])
        else:
            if eval(config['args']['save_comments']):
                shell("bwa mem -C -t "+str(threads)+" -R '@RG\\tID:"+params.sample+"\\tSM:"+params.sample+"' "+input.ref+" "+input.fq1+" "+input.fq2+" > "+output[0]+" 2> "+log[0])
            else:
                shell("bwa mem -t "+str(threads)+" -R '@RG\\tID:"+params.sample+"\\tSM:"+params.sample+"' "+input.ref+" "+input.fq1+" "+input.fq2+" > "+output[0]+" 2> "+log[0])

        

rule make_ref_te_bed:
    input:
        config['mcc']['formatted_ref_tes']
    
    threads: 1

    output:
        config['mcc']['ref_tes_bed']

    run:
        with open(input[0],"r") as i:
            with open(output[0], "w") as o:
                for line in i:
                    split_line = line.replace(";","\t").split("\t")
                    line = "\t".join([split_line[0], str(int(split_line[3])-1), split_line[4], split_line[9], ".", split_line[6]])
                    o.write(line+"\n")



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


rule retroseq:
    input:
        config['mcc']['ref_tes_bed'],
        config['mcc']['sam']
    
    output:
        config['args']['out']+"/retroseq/retroseq.log"
    
    run:
        shell("touch "+config['args']['out']+"/retroseq/retroseq.log")