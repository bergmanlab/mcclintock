
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
    input:
        ref = config['in']['reference'],
        consensus = config['in']['consensus'],

    params:
        coverage_fasta = config['in']['coverage_fasta']

    output:
        temp(config['mcc']['reference']),
        temp(config['mcc']['consensus']),
        temp(config['mcc']['coverage_fasta'])

    threads: 1
    
    run:
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+input.ref+" 80 > "+output[0])
        shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+input.consensus+" 80 > "+output[1])
        if params.coverage_fasta == "None":
            shell("touch "+output[2])
        else:
            shell("python "+config['args']['mcc_path']+"/scripts/fix_fasta_lines.py "+params.coverage_fasta+" 80 > "+output[2])


rule make_run_copies:
    output:
        temp(config['mcc']['locations']),
        temp(config['mcc']['taxonomy'])

    threads: 1

    run:
        if config['in']['locations'] == "None":
            shell("touch "+output[0])
        else:
            shell("cp "+config['in']['locations']+" "+output[0])
        
        if config['in']['taxonomy'] == "None":
            shell("touch "+output[1])
        else:
            shell("cp "+config['in']['taxonomy']+" "+output[1])

rule make_reference_te_files:
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
        config['mcc']['augmented_reference'],
        config['mcc']['popoolationTE_taxonomy'],
        config['mcc']['popoolationTE_gff']
    
    script:
        config['args']['mcc_path']+"/scripts/make_ref_te_files.py"

rule index_reference_genome:
    input:
        ref = config['mcc']['augmented_reference']
    
    threads: 1

    output:
        config['mcc']['augmented_reference']+".fai",
        config['mcc']['augmented_reference']+".bwt"
    
    run:
        shell("samtools faidx "+input.ref)
        shell("bwa index "+input.ref)

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

rule sam_to_bam:
    input:
        sam = config['mcc']['sam'],
        ref_idx = config['mcc']['augmented_reference']+".fai"
    
    threads: workflow.cores

    output:
        config['mcc']['bam'],
        config['summary']['flagstat'],
        temp(config['mcc']['mcc_files']+config['args']['run_id']+".tmp.bam")
    
    run:
        shell("samtools view -@ "+str(workflow.cores)+" -Sb -t "+input.ref_idx+" "+input.sam+" > "+output[2])
        shell("samtools sort -@ "+str(workflow.cores)+" "+output[2]+" "+output[0].replace(".bam", ""))
        shell("samtools index "+output[0])
        shell("samtools flagstat "+output[0]+" > "+output[1])
        

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

rule telocate_taxonomy:
    input:
        script = config['args']['mcc_path']+"/tools/te-locate/TE_hierarchy.pl",
        ref_gff = config['mcc']['formatted_ref_tes'],
        taxonomy = config['mcc']['formatted_taxonomy']
    
    threads: 1

    output:
        config['mcc']['telocate_te_gff']
    
    run:
        shell("perl "+input.script+" "+input.ref_gff+" "+input.taxonomy+" Alias")


rule median_insert_size:
    input:
        config['mcc']['sam']
    
    threads: 1

    output:
        config['summary']['median_insert_size']
    
    run:
        import statistics
        insert_sizes = []
        with open(input[0],"r") as sam:
            for line in sam:
                split_line = line.split("\t")
                if len(split_line) >= 8:
                    insert_size = int(split_line[8])
                    if insert_size > 0:
                        insert_sizes.append(insert_size)
        
        insert_sizes.sort()
        median = statistics.median(insert_sizes)
        with open(output[0],"w") as out:
            out.write("median_insert_size="+str(median)+"\n")

rule telocate_sam:
    input:
        config['mcc']['sam']
    
    threads: 1
    
    output:
        config['mcc']['telocate_sam']

    run:
        shell("sort -S "+config['args']['mem']+"G --temporary-directory="+config['args']['out']+"/tmp "+input[0]+" > "+output[0])

rule telocate_ref:
    input:
        config['mcc']['augmented_reference']
    
    threads: 1
    
    output:
        config['mcc']['telocate_ref_fasta']
    
    run:
        chromosomes = []
        with open(input[0],"r") as infa:
            with open(output[0],"w") as outfa:
                for line in infa:
                    outfa.write(line)
                    if ">" in line:
                        chromosomes.append(line)

        if len(chromosomes) < 5:
            diff = 5 - len(chromosomes)        
            with open(output[0],"a") as out: 
                for i in range(1,diff+1):
                    out.write(">fixforTElocate"+str(i)+"\n")
                    out.write("ACGT\n")

rule reference_2bit:
    input:
        config['mcc']['augmented_reference']
    
    output:
        config['mcc']['ref_2bit']
    
    run:
        shell("faToTwoBit "+input[0]+" "+output[0])

rule relocaTE_consensus:
    input:
        config['mcc']['formatted_consensus']
    
    threads: 1

    output:
        config['mcc']['relocaTE_consensus']

    run:
        with open(input[0],"r") as infa:
            with open(output[0], "w") as outfa:
                for line in infa:
                    if ">" in line:
                        line = line.replace("\n","")
                        line += " TSD=UNK\n"

                    outfa.write(line)
    
rule relocaTE_ref_gff:
    input:
        config['mcc']['formatted_ref_tes'],
        config['mcc']['formatted_taxonomy']

    threads: 1

    output:
        config['mcc']['relocaTE_ref_TEs']

    run:
        taxonomy = {}
        with open(input[1],"r") as tax:
            for line in tax:
                split_line = line.split("\t")
                taxonomy[split_line[0]] = split_line[1].replace("\n","")
        
        with open(input[0], "r") as ingff:
            with open(output[0],"w") as outgff:
                for line in ingff:
                    split_line = line.split("\t")
                    split_line[2] = taxonomy[split_line[2]]
                    line = "\t".join(split_line)
                    outgff.write(line)


rule popoolationTE_ref_fasta:
    input:
        config['mcc']['masked_fasta'],
        config['mcc']['formatted_consensus'],
        config['mcc']['ref_te_fasta']
    
    threads: 1

    output:
        config['mcc']['popoolationTE_ref_fasta']
    
    run:
        shell("cat "+input[0]+" "+input[1]+" "+input[2]+" > "+output[0])


rule coverage:
    input:
        config['mcc']['fq1'],
        config['mcc']['fq2'],
        config['mcc']['reference'],
        config['mcc']['consensus'],
        config['mcc']['coverage_fasta']
    
    threads: workflow.cores

    output:
        config['args']['out']+"/coverage/coverage.log"

    run:
        shell("touch "+config['args']['out']+"/coverage/coverage.log")    
    # script:
    #     "modules/coverage.py"

rule telocate:
    input:
        config['mcc']['telocate_te_gff'],
        config['summary']['median_insert_size'],
        config['mcc']['telocate_sam'],
        config['mcc']['telocate_ref_fasta'],
        config['args']['mcc_path']+"/install/tools/te-locate/TE_locate.pl"
    
    threads: 1

    output:
        config['args']['out']+"/te-locate/te-locate.log"
    
    run:
        shell("touch "+config['args']['out']+"/te-locate/te-locate.log")


rule retroseq:
    input:
        config['mcc']['formatted_consensus'],
        config['mcc']['bam'],
        config['mcc']['augmented_reference'],
        config['mcc']['ref_tes_bed'],
        config['mcc']['formatted_taxonomy']

    threads: 1
    
    output:
        config['args']['out']+"/retroseq/retroseq.log"
    
    run:
        shell("touch "+config['args']['out']+"/retroseq/retroseq.log")



rule TEMP:
    input:
        config['mcc']['bam'],
        config['mcc']['ref_2bit'],
        config['mcc']['formatted_consensus'],
        config['mcc']['ref_tes_bed'],
        config['mcc']['formatted_taxonomy'],
        config['summary']['median_insert_size'],
        config['mcc']['telocate_te_gff']

    threads: workflow.cores

    output:
        config['args']['out']+"/temp/temp.log"
    
    run:
        shell("touch "+output[0])



rule relocaTE:
    input:
        config['mcc']['relocaTE_consensus'],
        config['mcc']['relocaTE_ref_TEs'],
        config['mcc']['augmented_reference'],
        config['mcc']['fq1'],
        config['mcc']['fq2']

    threads: workflow.cores

    output:
        config['args']['out']+"/relocaTE/relocaTE.log"
    
    run:
        shell("touch "+output[0])


rule ngs_te_mapper:
    input:
        config['mcc']['formatted_consensus'],
        config['mcc']['augmented_reference'],
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: workflow.cores

    output:
        config['args']['out']+"/ngs_te_mapper/ngs_te_mapper.log"
    
    run:
        shell("touch "+output[0])
        

rule popoolationTE:
    input:
        config['mcc']['popoolationTE_ref_fasta'],
        config['mcc']['popoolationTE_taxonomy'],
        config['mcc']['popoolationTE_gff'],
        config['mcc']['fq1'],
        config['mcc']['fq2']
    
    threads: workflow.cores

    output:
        config['args']['out']+"/popoolationTE/popoolationTE.log"

    run:
        shell("touch "+output[0])

