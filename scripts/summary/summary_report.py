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

class Prediction:
    def __init__(self):
        self.family = ""
        self.all = []
        self.reference = []
        self.nonreference = []

class MethodPrediction:
    def __init__(self):
        self.method = ""
        self.family = ""
        self.all = []
        self.reference = []
        self.nonreference = []
        self.insertions = []

class Insertion:
    def __init__(self):
        self.chrom = ""
        self.family = ""
        self.type = ""
        self.start = -1
        self.end = -1
        self.strand = "."


def main():
    out_files = snakemake.input.out_files
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2

    commit = snakemake.params.commit
    ref = snakemake.params.ref
    consensus = snakemake.params.consensus
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
    chromosomes = snakemake.params.chromosomes.split(",")
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
    make_families_page(env, consensus, methods, out_file_map, out_dir)
    make_family_pages(env, consensus, methods, out_file_map, chromosomes, out_dir)

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
    mccutils.mkdir(out_dir+"/html/")
    mccutils.mkdir(out_dir+"/css/")
    for css in os.listdir(css_dir):
        mccutils.run_command(["cp", css_dir+"/"+css, out_dir+"/css/"])
    
    mccutils.mkdir(out_dir+"/js/")
    for js in os.listdir(js_dir):
        mccutils.run_command(["cp", js_dir+"/"+js, out_dir+"/js/"])

def make_data_copies(methods, results_dir, out_dir):
    mccutils.mkdir(out_dir+"/data/")
    if "trimgalore" in methods:
        mccutils.mkdir(out_dir+"/data/trimgalore/")
        for f in os.listdir(results_dir+"/trimgalore"):
            if ".zip" not in f:
                mccutils.run_command(["cp", "-r", results_dir+"/trimgalore/"+f, out_dir+"/data/trimgalore/"])
    
    if "coverage" in methods:
        mccutils.mkdir(out_dir+"/data/coverage/")
        for f in os.listdir(results_dir+"/coverage/"):
            if not os.path.isdir(results_dir+"/coverage/"+f):
                mccutils.run_command(["cp", results_dir+"/coverage/"+f, out_dir+"/data/coverage/"])
        for f in os.listdir(results_dir+"/coverage/te-depth-files/"):
            mccutils.run_command(["cp", results_dir+"/coverage/te-depth-files/"+f, out_dir+"/data/coverage/"])
        for f in os.listdir(out_dir+"/data/coverage/"):
            o = f.replace(".csv",".txt")
            o = o.replace(".cov",".txt")
            mccutils.run_command(["mv", out_dir+"/data/coverage/"+f, out_dir+"/data/coverage/"+o])

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


def make_families_page(jinja_env, consensus, methods, out_file_map, out_dir):
    template = jinja_env.get_template('families.html')

    families = []
    with open(consensus,"r") as fa:
        for line in fa:
            if line[0] == ">":
                family = line.replace(">","")
                family = family.replace("\n","")
                families.append(family)
    
    prediction_methods = []
    for method in methods:
        if method != "coverage" and method != "trimgalore":
            prediction_methods.append(method)

    prediction_list = []
    for family in families:
        prediction = count_predictions(prediction_methods,out_file_map, family)
        prediction_list.append(prediction)

    mccutils.mkdir(out_dir+"/data/families/")
    prediction_summary_file = out_dir+"/data/families/family_prediction_summary.txt"
    write_prediction_file(prediction_list, prediction_methods, prediction_summary_file)

    rendered_lines = template.render(
        methods=prediction_methods,
        families=families,
        predictions=prediction_list
    )

    out_file = out_dir+"/html/families.html"
    with open(out_file,"w") as out:
        for line in rendered_lines:
            out.write(line)


def make_family_pages(jinja_env, consensus, methods, out_file_map, chromosomes, out_dir):
    families = []
    with open(consensus,"r") as fa:
        for line in fa:
            if line[0] == ">":
                family = line.replace(">","")
                family = family.replace("\n","")
                families.append(family)

    prediction_methods = []
    for method in methods:
        if method != "coverage" and method != "trimgalore":
            prediction_methods.append(method)


    depth = {}
    with open(out_dir+"/data/coverage/te_depth.txt","r") as depth_file:
        for i,line in enumerate(depth_file):
            if i > 0:
                split_line = line.split(",")
                family = split_line[0]
                depth[family] = [split_line[1], split_line[2]]

    for family in families:
        all_cov = None
        all_pos = None
        uniq_cov = None
        uniq_pos = None

        template = jinja_env.get_template('family.html')

        if "coverage" in methods:
            all_pos = []
            all_cov = []
            with open(out_dir+"data/coverage/"+family+".allQ.normalized.txt","r") as data:
                for line in data:
                    line = line.replace("\n","")
                    split_line = line.split("\t")
                    all_pos.append(split_line[1])
                    all_cov.append(split_line[2])
            
            uniq_pos = []
            uniq_cov = []
            with open(out_dir+"data/coverage/"+family+".highQ.normalized.txt","r") as data:
                for line in data:
                    line = line.replace("\n","")
                    split_line = line.split("\t")
                    uniq_pos.append(split_line[1])
                    uniq_cov.append(split_line[2])

        prediction = count_predictions(prediction_methods, out_file_map, family)
        family_prediction_summary_file = out_dir+"/data/families/"+family+"_prediction_summary.txt"
        write_prediction_file([prediction], prediction_methods, family_prediction_summary_file)

        method_predictions = []
        for method in prediction_methods:
            method_prediction = count_predictions_chrom(method, out_file_map, family, chromosomes)
            method_prediction.insertions = get_family_predictions(out_file_map[method], family)
            family_predictions_file = out_dir+"/data/families/"+family+"_"+method+"_predictions.txt"
            with open(family_predictions_file,"w") as predictions_file:
                for insertion in method_prediction.insertions:
                    line = ",".join([insertion.chrom, insertion.family, insertion.type, str(insertion.start), str(insertion.end), insertion.strand])
                    predictions_file.write(line+"\n")
            method_predictions.append(method_prediction)

        rendered_lines = template.render(
            methods=prediction_methods,
            family=family,
            all_coverage=all_cov,
            all_positions=all_pos,
            uniq_coverage=uniq_cov,
            uniq_positions=uniq_pos,
            norm_depth=depth[family][0],
            uniq_depth=depth[family][1],
            prediction_summary=prediction,
            chromosomes=chromosomes,
            method_results=method_predictions
        )
        
        out_file = out_dir+"/html/"+family+".html"
        with open(out_file,"w") as out:
            for line in rendered_lines:
                out.write(line)


def count_predictions(methods, out_file_map, family):
        prediction = Prediction()
        prediction.family = family

        for method in methods:
            all_count = 0
            reference_count = 0
            nonreference_count = 0
            prediction_file = out_file_map[method]
            with open(prediction_file,"r") as predictions:
                for line in predictions:
                    split_line = line.split("\t")
                    if len(split_line) > 3:
                        info = split_line[3]
                        if "_reference_" in info:
                            split_info = info.split("_reference_")
                        else:
                            split_info = info.split("_non-reference_")
                        if split_info[0] == family:
                            all_count += 1
                            if "non-reference" in info:
                                nonreference_count += 1
                            else:
                                reference_count += 1
            
            prediction.all.append(all_count)
            prediction.reference.append(reference_count)
            prediction.nonreference.append(nonreference_count)
        
        return prediction

def count_predictions_chrom(method, out_file_map, family, chromosomes):
        prediction = MethodPrediction()
        prediction.method = method
        prediction.family = family

        for chromosome in chromosomes:
            all_count = 0
            reference_count = 0
            nonreference_count = 0
            prediction_file = out_file_map[method]
            with open(prediction_file,"r") as predictions:
                for line in predictions:
                    split_line = line.split("\t")
                    if len(split_line) > 3:
                        info = split_line[3]
                        if "_reference_" in info:
                            split_info = info.split("_reference_")
                        else:
                            split_info = info.split("_non-reference_")
                        if split_info[0] == family and split_line[0] == chromosome:
                            all_count += 1
                            if "non-reference" in info:
                                nonreference_count += 1
                            else:
                                reference_count += 1
            
            prediction.all.append(all_count)
            prediction.reference.append(reference_count)
            prediction.nonreference.append(nonreference_count)
        
        return prediction

def get_family_predictions(bed, family):
    predictions = []
    with open(bed,"r") as infile:
        for line in infile:
            line = line.replace("\n","")
            split_line = line.split("\t")
            if len(split_line) > 3:
                info = split_line[3]
                insert_type = ""
                if "_reference_" in info:
                    split_info = info.split("_reference_")
                    insert_type = "Reference"
                else:
                    split_info = info.split("_non-reference_")
                    insert_type = "Non-Reference"
                if split_info[0] == family:
                    insertion = Insertion()
                    insertion.chrom = split_line[0]
                    insertion.family = family
                    insertion.start = int(split_line[1])
                    insertion.end = int(split_line[2])
                    insertion.strand = split_line[5]
                    insertion.type = insert_type
                    predictions.append(insertion)
    
    return predictions



def write_prediction_file(prediction_list, methods, out_file):
    with open(out_file,"w") as out:
        header = ["TE_Family","Type"] + methods
        header = ",".join(header)
        out.write(header+"\n")
        for prediction in prediction_list:
            line = [prediction.family,"all"]
            for val in prediction.all:
                line.append(str(val))
            line = ",".join(line)
            out.write(line+"\n")

            line = [prediction.family,"reference"]
            for val in prediction.reference:
                line.append(str(val))
            line = ",".join(line)
            out.write(line+"\n")

            line = [prediction.family,"non-reference"]
            for val in prediction.nonreference:
                line.append(str(val))
            line = ",".join(line)
            out.write(line+"\n")

if __name__ == "__main__":                
    main()

