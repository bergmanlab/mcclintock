import os
import sys
import subprocess
from Bio import SeqIO
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.config as config
from datetime import datetime
from jinja2 import Environment, FileSystemLoader

templateLoader = FileSystemLoader(searchpath=snakemake.config['args']['mcc_path']+"/templates/html/")
env = Environment(loader=templateLoader)


def main():
    out_files = snakemake.input.out_files
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2

    commit = snakemake.params.commit
    ref = snakemake.params.ref
    taxonomy = snakemake.params.taxonomy
    bam = snakemake.params.bam
    flagstat = snakemake.params.flagstat
    median_insert_size = snakemake.params.median_insert_size
    methods = snakemake.params.methods
    results_dir = snakemake.params.results_dir
    sample_name = snakemake.params.sample_name
    command = snakemake.params.command
    execution_dir = snakemake.params.execution_dir
    start_time = snakemake.params.time
    raw_fq2 = snakemake.params.raw_fq2
    out_dir = snakemake.params.out_dir

    paired = True
    if raw_fq2 == "None":
        paired = False

    out_file_map = {}
    for method in config.OUT_PATHS.keys():
        out_file_map[method] = config.OUT_PATHS[method]
        out_file_map[method] = out_file_map[method].replace("{{results}}", results_dir)
        out_file_map[method] = out_file_map[method].replace("{{samplename}}", sample_name)

    if os.path.exists(taxonomy):
        make_te_csv(methods, out_file_map, taxonomy, out_dir+"te_summary.csv")

    mapping_info,end_time = make_run_summary(out_file_map, commit, methods, fq1, fq2, ref, bam, flagstat, median_insert_size, command, execution_dir, start_time, out_dir, snakemake.output.summary_report, paired=paired)
    make_local_css_js_copies(snakemake.config['args']['mcc_path']+"/templates/css/", snakemake.config['args']['mcc_path']+"/templates/js/", snakemake.params.out_dir)
    make_data_copies(methods, snakemake.params.results_dir, snakemake.params.out_dir)
    make_summary_page(env, methods, sample_name, commit, start_time, end_time, out_dir, execution_dir, command, snakemake.input.fq1, snakemake.input.fq2, mapping_info, out_file_map, paired, snakemake.output.html_summary_report)
    


def get_te_counts(bed):
    all_te = 0
    ref_te = 0
    nonref_te = 0

    with open(bed,"r") as b:
        for line in b:
            split_line = line.split("\t")
            if len(split_line) == 6:
                all_te += 1
                if "non-reference" in split_line[3]:
                    nonref_te += 1
                else:
                    ref_te += 1


    return all_te, ref_te, nonref_te


def pad(string, total_len, symbol=" ", front=False):
    if not front:
        string = string+(symbol*(total_len - len(string)))
    else:
        string = (symbol*(total_len - len(string)))+string

    return string
    


def make_te_csv(methods, out_file_map, taxonomy, out_csv):
    te_counts = {}
    method_counts = {}
    header = ["TE-Family"]

    # get all TE family names
    te_names = []
    with open(taxonomy,"r") as taxon:
        for line in taxon:
            te_name = line.split("\t")[1].replace("\n","")
            if te_name not in te_names:
                te_names.append(te_name)
    
    # setup TE family count dicts and header
    for x,te_name in enumerate(te_names):
        te_counts[te_name] = {}
        for method in config.ALL_METHODS:
            if "nonredundant.bed" in out_file_map[method]:
                if method in methods:
                    te_counts[te_name][method] = [0,0]
                else:
                    te_counts[te_name][method] = ["NA","NA"]
                
                if x < 1:
                    header += [method+"_all", method+"_ref", method+"_nonref"]
    
    
    # count TE family predictions for each method
    for method in methods:
        bed = out_file_map[method]
        if "nonredundant.bed" in bed:
            with open(bed, "r") as b:
                for line in b:
                    split_line = line.split("\t")
                    if len(split_line) == 6:
                        name = split_line[3]
                        te_name = ""
                        split_name = name.split("_")
                        for x,part in enumerate(split_name):
                            if part == "reference" or part == "non-reference":
                                te_name = "_".join(split_name[:x])
                        
                        if te_name == "":
                            print("ERROR: couldn't find TE name in:", name)
                    
                        else:
                            if te_name not in te_counts.keys():
                                te_counts[te_name] = {method:[0,0]}
                            elif method not in te_counts[te_name].keys():
                                te_counts[te_name][method] = [0,0]
                                
                            if "non-reference" in name:
                                te_counts[te_name][method][1] += 1
                            else:
                                te_counts[te_name][method][0] += 1

    # create te family count csv
    with open(out_csv,"w") as out:
        out.write(",".join(header)+"\n")
        for te in te_names:
            line = [te.replace(",","_")]
            for method in config.ALL_METHODS:
                if "nonredundant.bed" in out_file_map[method]:
                    if method in methods:
                        line += [str(te_counts[te][method][0]+te_counts[te][method][1]), str(te_counts[te][method][0]), str(te_counts[te][method][1])]
                    else:
                        line +=["NA","NA","NA"]
            
            out.write(",".join(line)+"\n")
    

def make_run_summary(out_file_map, commit, methods, fq1, fq2, ref, bam, flagstat, median_insert_size, command, execution_dir, start_time, out_dir, out_file, paired=False):
    out_lines = ["\n"]
    out_lines.append(("-"*34)+"\n")
    out_lines.append("MCCLINTOCK SUMMARY REPORT\n")
    out_lines.append(("-"*34) + "\n")
    split_command = command.split(" ")
    for x,split in enumerate(split_command):
        if split[0] == "-":
            split_command[x] = split_command[x].replace("-", " \\\n\t-",1)
    command = " ".join(split_command)
    out_lines.append("McClintock Version: "+commit+"\n\n")
    out_lines.append("Command:\n"+command+"\n")
    out_lines.append("\nrun from directory: "+execution_dir+"\n")
    out_lines.append(pad("Started:",12)+start_time+"\n")
    out_lines.append(pad("Completed:",12)+"{{END_TIME}}"+"\n\n")

    if os.path.exists(bam) and os.path.exists(flagstat) and os.path.exists(median_insert_size) and os.path.exists(ref):
        mapping_info = {}
        out_lines.append(("-"*34)+"\n")
        out_lines.append("MAPPED READ INFORMATION\n")
        out_lines.append(("-"*34) + "\n")

        fq1_read_len = str(mccutils.estimate_read_length(fq1))
        mapping_info["read1_length"] = fq1_read_len
        out_lines.append(pad("read1 sequence length:",24) + fq1_read_len+"\n")
        if paired:
            fq2_read_len = str(mccutils.estimate_read_length(fq2))
            mapping_info["read2_length"] = fq2_read_len
            out_lines.append(pad("read2 sequence length:",24) + fq2_read_len+"\n")
        
        with open(flagstat,"r") as stat:
            for line in stat:
                if "read1" in line:
                    reads = line.split("+")[0].replace(" ","")
                    mapping_info['read1_reads'] = str(reads)
                    out_lines.append(pad("read1 reads:",24) + str(reads) + "\n")
                
                if paired and "read2" in line:
                    reads = line.split("+")[0].replace(" ","")
                    mapping_info['read2_reads'] = str(reads)
                    out_lines.append(pad("read2 reads:",24) + str(reads) + "\n")

            
        with open(median_insert_size,"r") as median_insert:
            for line in median_insert:
                line = line.split("=")[1].replace("\n","")
                insert_size = str(int(float(line)))
                mapping_info["median_insert_size"] = insert_size
                out_lines.append(pad("median insert size:",24) + insert_size + "\n")
        
        avg_genome_cov = str(get_avg_coverage(ref, bam, out_dir))
        mapping_info['avg_genome_cov'] = avg_genome_cov
        out_lines.append(pad("avg genome coverage:",24) + avg_genome_cov + "\n")
        out_lines.append(("-"*34)+"\n")


    len_longest_name = 0
    for method in config.ALL_METHODS:
        if len(method) > len_longest_name:
            len_longest_name = len(method)


    width1 = len_longest_name+2
    width2 = 10
    width3 = 13
    width4 = 14

    out_lines.append("\n")
    out_lines.append("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")
    out_lines.append(pad("METHOD", width1) + pad("ALL",width2) + pad("REFERENCE",width3) + pad("NON-REFERENCE", width4) + "\n")
    out_lines.append("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")
    for method in config.ALL_METHODS:
        if "nonredundant.bed" in out_file_map[method]:
            if method in methods:
                all_te, ref_te, nonref_te = get_te_counts(out_file_map[method])
                out_lines.append(pad(method, width1) + pad(str(all_te), width2) + pad(str(ref_te), width3) + pad(str(nonref_te), width4) + "\n")
            else:
                out_lines.append(pad(method, width1) + pad("NA", width2) + pad("NA", width3) + pad("NA", width4) + "\n")

    out_lines.append("-"*(width1) + "-"*width2 + "-"*width3 + "-"*width4 + "\n")

        
    with open(out_file,"w") as out:
        for line in out_lines:
            if "{{END_TIME}}" in line:
                now = datetime.now()
                completed = now.strftime("%Y-%m-%d %H:%M:%S")
                line = line.replace("{{END_TIME}}",completed)
            print(line, end="")
            out.write(line)
    
    return mapping_info,completed


def get_avg_coverage(ref, bam, out):
    chrom = []
    fasta_records = SeqIO.parse(ref,"fasta")
    for record in fasta_records:
        chrom.append(str(record.id))
    
    tmp = out+"/tmp"
    command = ['samtools','depth',bam]
    mccutils.run_command_stdout(command, tmp)

    cov_total = 0
    pos = 0 
    with open(tmp,"r") as depth:
        for line in depth:
            split_line = line.split("\t")
            if split_line[0] in chrom:
                pos += 1
                cov_total += int(split_line[2])
    
    mccutils.remove(tmp)
    return round(cov_total/pos, 3)

def make_local_css_js_copies(css_dir, js_dir, out_dir):
    mccutils.mkdir(out_dir+"/css/")
    for css in os.listdir(css_dir):
        mccutils.run_command(["cp", css_dir+"/"+css, out_dir+"/css/"])
    
    mccutils.mkdir(out_dir+"/js/")
    for js in os.listdir(js_dir):
        mccutils.run_command(["cp", js_dir+"/"+js, out_dir+"/js/"])

def make_data_copies(methods, results_dir, out_dir):
    if "trimgalore" in methods:
        mccutils.mkdir(out_dir+"/data/")
        mccutils.mkdir(out_dir+"/data/trimgalore/")
        for f in os.listdir(results_dir+"/trimgalore"):
            if ".zip" not in f:
                mccutils.run_command(["cp", "-r", results_dir+"/trimgalore/"+f, out_dir+"/data/trimgalore/"])

def read_trimgalore_results(fastq, trimgalore_dir):
    results = ["","","","","",""]
    trimgalore_file = ""
    fq_base_name = fastq.split("/")[-1]
    for f in os.listdir(trimgalore_dir):
        if ".txt" in f and fq_base_name in f:
            trimgalore_file = f
            with open(trimgalore_dir+"/"+trimgalore_file,"r") as res:
                for line in res:
                    if "Total reads processed:" in line:
                        line = line.split(":")[1]
                        line = line.replace(" ","")
                        line = line.replace(",","")
                        results[0] = line
                    
                    elif "Reads with adapters:" in line:
                        line = line.split(":")[1]
                        line = line.replace(" ","")
                        line = line.replace("(", " (")
                        line = line.replace(",","")
                        results[1] = line
                    
                    elif "Reads written (passing filters):" in line:
                        line = line.split(":")[1]
                        line = line.replace(" ","")
                        line = line.replace("(", " (")
                        line = line.replace(",","")
                        results[2] = line
                    
                    elif "Total basepairs processed:" in line:
                        line = line.split(":")[1]
                        line = line.replace(" ","")
                        line = line.replace("b", " b")
                        line = line.replace(",","")
                        results[3] = line
                    
                    elif "Quality-trimmed:" in line:
                        line = line.split(":")[1]
                        line = line.replace(" ","")
                        line = line.replace("b", " b")
                        line = line.replace("(", " (")
                        line = line.replace(",","")
                        results[4] = line
                    
                    elif "Total written (filtered):" in line:
                        line = line.split(":")[1]
                        line = line.replace(" ","")
                        line = line.replace("b", " b")
                        line = line.replace("(", " (")
                        line = line.replace(",","")
                        results[5] = line
    
    return results, trimgalore_file


def make_summary_page(jinja_env, methods, sample_name, commit, start_time, end_time, out_dir, execution_dir, command, fq1, fq2, mapping_info, out_file_map, paired, out_file):
    template = jinja_env.get_template('summary.html')

    # split command into separate lines
    split_command = command.split(" ")
    for x,split in enumerate(split_command):
        if split[0] == "-":
            split_command[x] = split_command[x].replace("-", " \{{break}}-",1)
    command = " ".join(split_command)

    split_command = command.split("{{break}}")

    # get trimgalore results
    if "trimgalore" in methods:
        fq1_trimgalore_results, fq1_trimgalore_file = read_trimgalore_results(fq1, out_dir+"/data/trimgalore/")
        fq1_trimgalore_file = "data/trimgalore/"+fq1_trimgalore_file

        if paired:
            fq2_trimgalore_results, fq2_trimgalore_file = read_trimgalore_results(fq2, out_dir+"/data/trimgalore/")
            fq2_trimgalore_file = "data/trimgalore/"+fq2_trimgalore_file
        else:
            fq2_trimgalore_results = None
    else:
        fq1_trimgalore_results = None
        fq2_trimgalore_results = None


    # get mapping info
    read1_seq_len = mapping_info['read1_length']
    read2_seq_len = None
    read1_reads = mapping_info['read1_reads']
    read2_reads = None
    median_insert_size = mapping_info['median_insert_size']
    avg_genome_cov = mapping_info['avg_genome_cov']
    if paired:
        read2_seq_len = mapping_info['read2_length']
        read2_reads = mapping_info['read2_reads']
        

    # get prediction information
    prediction_methods = []
    reference_counts = []
    nonreference_counts = []
    for method in methods:
        if method != "trimgalore" and method != "coverage":
            prediction_methods.append(method)
            results = out_file_map[method]
            nonreference = 0
            reference = 0
            with open(results,"r") as bed:
                for line in bed:
                    split_line = line.split("\t")
                    if len(split_line) > 3:
                        if "_reference_" in split_line[3]:
                            reference += 1
                        elif "_non-reference_" in split_line[3]:
                            nonreference += 1
            reference_counts.append(reference)
            nonreference_counts.append(nonreference)
            

    with open(out_dir+"/data/run/te_prediction_summary.txt","w") as out:
        out.write("method,reference,non-reference\n")
        for i,method in enumerate(prediction_methods):
            out_line = ",".join([method, str(reference_counts[i]), str(nonreference_counts[i])])
            out.write(out_line+"\n")


    rendered_lines = template.render(
        paired=paired,
        sample=sample_name,
        methods=methods,
        commit=commit,
        run_start=start_time,
        run_end=end_time,
        run_path=execution_dir,
        command=split_command,
        trimgalore1=fq1_trimgalore_results,
        trimgalore2=fq2_trimgalore_results,
        trimgalore1_file=fq1_trimgalore_file,
        trimgalore2_file=fq2_trimgalore_file,
        r1_seq_len = read1_seq_len,
        r2_seq_len = read2_seq_len,
        r1_reads = read1_reads,
        r2_reads = read2_reads,
        median_insert_size = median_insert_size,
        avg_genome_cov = avg_genome_cov,
        prediction_methods=prediction_methods,
        reference_counts=reference_counts,
        nonreference_counts=nonreference_counts
    )

    with open(out_file,"w") as out:
        for line in rendered_lines:
            out.write(line)

if __name__ == "__main__":                
    main()

