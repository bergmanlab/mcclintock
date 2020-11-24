import argparse
import os
import sys
import traceback
import json
import random
import subprocess
import statistics
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict

class Insertion:
    def __init__(self):
        self.chromosome = ""
        self.family = ""
        self.start = -1
        self.end = -1
        self.strand = "?"
        self.reference = False

def main():
    args = parse_args()
    if args.mode == "run":
        for x in range(args.start,args.end+1):
            # forward
            consensus_seqs = get_seqs(args.consensus)
            reference_seqs = get_seqs(args.reference)
            
            modified_reference = args.out+"/data/forward/"+str(args.runid)+str(x)+".modref.fasta"
            if not os.path.exists(modified_reference):
                modified_reference = add_synthetic_insertion(reference_seqs, consensus_seqs, args.config, x, args.out, run_id=args.runid, seed=args.seed)
            
            fastq1 = modified_reference.replace(".fasta", "_1.fastq")
            fastq2 = modified_reference.replace(".fasta", "_2.fastq")
            if not os.path.exists(fastq1) or not os.path.exists(fastq2):
                num_pairs = calculate_num_pairs(modified_reference)
                fastq1, fastq2 = create_synthetic_reads(modified_reference, num_pairs, x, args.out, run_id=args.runid, seed=args.seed)

            if not os.path.exists(args.out+"/results/forward/run"+args.runid+"_"+str(x)+"/results/summary/summary_report.txt"):
                run_mcclintock(fastq1, fastq2, args.reference, args.consensus, args.locations, args.taxonomy, x, args.proc, args.out, args.config, run_id=args.runid)

            # reverse
            consensus_seqs = get_seqs(args.consensus)
            reference_seqs = get_seqs(args.reference)

            modified_reference = args.out+"/data/reverse/"+str(args.runid)+str(x)+".modref.fasta"
            if not os.path.exists(modified_reference):
                modified_reference = add_synthetic_insertion(reference_seqs, consensus_seqs, args.config, x, args.out, run_id=args.runid, seed=args.seed, reverse=True)
            
            fastq1 = modified_reference.replace(".fasta", "_1.fastq")
            fastq2 = modified_reference.replace(".fasta", "_2.fastq")
            if not os.path.exists(fastq1) or not os.path.exists(fastq2):
                num_pairs = calculate_num_pairs(modified_reference)
                fastq1, fastq2 = create_synthetic_reads(modified_reference, num_pairs, x, args.out, run_id=args.runid, seed=args.seed)

            if not os.path.exists(args.out+"/results/reverse/run"+args.runid+"_"+str(x)+"/results/summary/summary_report.txt"):
                run_mcclintock(fastq1, fastq2, args.reference, args.consensus, args.locations, args.taxonomy, x, args.proc, args.out, args.config, run_id=args.runid, reverse=True)

    elif args.mode == "analysis":
        if not os.path.exists(args.out+"/summary/"):
            os.mkdir(args.out+"/summary/")

        actual_insertions = get_actual_insertions(args.out+"/data/forward/")
        predicted_insertions, methods = get_predicted_insertions(args.out+"/results/forward/")
        make_out_table(actual_insertions, predicted_insertions, methods, args.out+"/summary/forward_summary.csv")

        actual_insertions = get_actual_insertions(args.out+"/data/reverse/")
        predicted_insertions, methods = get_predicted_insertions(args.out+"/results/reverse/")
        make_out_table(actual_insertions, predicted_insertions, methods, args.out+"/summary/reverse_summary.csv")


def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock Simulation', description="Script to run synthetic insertion simulations to evaluate mcclintock component methods")

    ## required ##
    parser.add_argument("--mode", type=str, help="which mode this script should be run", required=True)
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required='run' in sys.argv)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required='run' in sys.argv)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required='run' in sys.argv)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required='run' in sys.argv)
    parser.add_argument("-j","--config", type=str, help="A json config file containing information on TE family TSD size and target sites", required='run' in sys.argv)
    ## optional ##
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required='run' not in sys.argv)
    parser.add_argument("--start", type=int, help="The number of replicates to run. [default = 1]", required=False)
    parser.add_argument("--end", type=int, help="The number of replicates to run. [default = 300]", required=False)
    parser.add_argument("--seed", type=str, help="a seed to the random number generator so runs can be replicated", required=False)
    parser.add_argument("--runid", type=str, help="a string to prepend to output files so that multiple runs can be run at the same time without file name clashes", required=False)

    args = parser.parse_args()

    # check --mode
    valid_modes = ["run", "analysis"]
    if args.mode not in valid_modes:
        print(args.mode, "not a valid mode ( valid modes:",",".join(valid_modes),")", file=sys.stderr)

    if args.mode == "run":
        #check -r
        args.reference = get_abs_path(args.reference)
        #check -c
        args.consensus = get_abs_path(args.consensus)
        # check -g
        args.locations = get_abs_path(args.locations)
        # check -t
        args.taxonomy = get_abs_path(args.taxonomy)
        # check -j
        args.config = get_abs_path(args.config)
        with open(args.config, "r") as j:
            config = json.load(j, object_pairs_hook = OrderedDict)
        args.config = config

        #check -p
        if args.proc is None:
            args.proc = 1

        #check -o
        if args.out is None:
            args.out = os.path.abspath(".")
        else:
            args.out = os.path.abspath(args.out)

            if not os.path.exists(args.out):
                try:
                    os.mkdir(args.out)
                except Exception as e:
                    track = traceback.format_exc()
                    print(track, file=sys.stderr)
                    print("cannot create output directory: ",args.out,"exiting...", file=sys.stderr)
                    sys.exit(1)

        # check --start
        if args.start is None:
            args.start = 1

        # check --end
        if args.end is None:
            args.end = 1
        
        if args.runid is None:
            args.runid = ""

    return args

def run_command(cmd_list, log=None):
    msg = ""
    if log is None:
        try:
            # print(" ".join(cmd_list))
            subprocess.check_call(cmd_list)
        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            sys.stderr.write(msg)
            sys.exit(1)
    
    else:
        try:
            out = open(log,"a")
            out.write(" ".join(cmd_list)+"\n")
            subprocess.check_call(cmd_list, stdout=out, stderr=out)
            out.close()

        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            writelog(log, msg)
            sys.stderr.write(msg)
            sys.exit(1)

def run_command_stdout(cmd_list, out_file, log=None):
    msg = ""
    if log is None:
        try:
            # print(" ".join(cmd_list)+" > "+out_file)
            out = open(out_file,"w")
            subprocess.check_call(cmd_list, stdout=out)
            out.close()
        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            sys.stderr.write(msg)
            sys.exit(1)
    
    else:
        try:
            out_log = open(log,"a")
            out_log.write(" ".join(cmd_list)+" > "+out_file+"\n")
            out = open(out_file,"w")
            subprocess.check_call(cmd_list, stdout=out, stderr=out_log)
            out.close()
            out_log.close()

        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            writelog(log, msg)
            sys.stderr.write(msg)
            sys.exit(1)

def writelog(log, msg):
    if log is not None:
        with open(log, "a") as out:
            out.write(msg)

def get_abs_path(in_file, log=None):
    if os.path.isfile(in_file):
        return os.path.abspath(in_file)
    else:
        msg = " ".join(["Cannot find file:",in_file,"exiting....\n"])
        sys.stderr.write(msg)
        writelog(log, msg)
        sys.exit(1)

def get_seqs(fasta):
    seq_dir = []
    records = SeqIO.parse(fasta, "fasta")
    for record in records:
        seq_dir.append([str(record.id), str(record.seq)])
    
    return seq_dir

def fix_fasta_lines(fasta, length):
    lines = []
    fasta_records = SeqIO.parse(fasta,"fasta")
    for record in fasta_records:
        # print(">"+record.id)
        header = ">"+str(record.id)
        lines.append(header)
        seq = str(record.seq)
        x = 0
        while(x+length < len(seq)):
            # print(seq[x:x+length])
            lines.append(seq[x:x+length])
            x += length

        remainder = (len(seq)) - x
        # print(seq[x:x+remainder])
        lines.append(seq[x:x+remainder])
    
    return lines

def add_synthetic_insertion(reference_seqs, consensus_seqs, config, rep, out, run_id="", seed=None, reverse=False):
    if not os.path.exists(out+"/data"):
        os.mkdir(out+"/data")
    
    if not reverse:
        outdir = out+"/data/forward"
    else:
        outdir = out+"/data/reverse"
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if seed is not None:
        random.seed(seed+"add_synthetic_"+str(reverse)+"_insertion"+str(rep))
    else:
        random.seed(str(datetime.now())+"add_synthetic_"+str(reverse)+"_insertion"+str(rep)+str(reverse))


    # get TE family to add
    te_families = list(config["families"].keys())
    family_to_add = te_families[random.randint(0, len(te_families)-1)]
    family_seq = ""
    for seq in consensus_seqs:
        if seq[0] == family_to_add:
            if not reverse:
                family_seq = seq[1]
            else:
                tmp = Seq(seq[1])
                tmp = tmp.reverse_complement()
                family_seq = str(tmp)
    target_file = config['families'][family_to_add]["targets"]
    valid_chroms = []
    with open(target_file,"r") as t:
        for line in t:
            chrom = line.split("\t")[0]
            valid_chroms.append(chrom)
    
    # get target site to add TE
    targets = []
    with open(target_file,"r") as t:
        for line in t:
            split_line = line.split("\t")
            targets.append(split_line)

    if len(targets) == 0:
        sys.exit("no targets in target bed file...exiting...")
    
    target = targets[random.randint(0, len(targets)-1)]
    target_site = random.randint(int(target[1]), int(target[2])-1)
    target_chrom = target[0]
    duplication_size = config['families'][family_to_add]["TSD"]


    chrom_idx_to_change = None
    seq_to_change = None

    # get reference contig to modify
    for x in range(0,len(reference_seqs)):
        if reference_seqs[x][0] == target_chrom:
            chrom_idx_to_change = x
            seq_to_change = reference_seqs[x][1]

    if chrom_idx_to_change is None:
        sys.exit("chrom: "+target_chrom+" not in reference")

    # add TE in reference contig
    seq_start = seq_to_change[0:target_site+duplication_size]
    seq_end = seq_to_change[target_site:]
    modified_seq = seq_start + family_seq + seq_end
    reference_seqs[chrom_idx_to_change][1] = modified_seq


    outfasta = outdir+"/"+run_id+str(rep)+".modref.fasta"
    outbed = outdir+"/"+run_id+str(rep)+".modref.bed"

    
    with open(outfasta+".tmp","w") as outfa:
        for seq in reference_seqs:
            outfa.write(">"+seq[0]+"\n")
            outfa.write(seq[1]+"\n")
    
    # create fasta with fixed line lengths
    fasta_lines = fix_fasta_lines(outfasta+".tmp", 80)
    with open(outfasta,"w") as outfa:
        for line in fasta_lines:
            outfa.write(line+"\n")

    os.remove(outfasta+".tmp")
    
    # create bed that marks TSD
    with open(outbed,"w") as bed:
        bed.write(target_chrom+"\t"+str(target_site)+"\t"+str(target_site+duplication_size)+"\t"+family_to_add+"\n")

    print(target_chrom+"\t"+str(target_site)+"\t"+str(target_site+duplication_size)+"\t"+family_to_add+"\n")

    return outfasta


def calculate_num_pairs(fasta):
    command = ["samtools","faidx", fasta]
    run_command(command)

    total_length = 0
    with open(fasta+".fai", "r") as faidx:
        for line in faidx:
            split_line = line.split("\t")
            contig_length = int(split_line[1])
            total_length += contig_length
    
    num_pairs = (total_length * 100)/202
    return num_pairs

def create_synthetic_reads(reference, num_pairs, rep, out, run_id="", seed=None):
    if seed is not None:
        random.seed(seed+"create_synthetic_reads"+str(rep))
    else:
        random.seed(str(datetime.now())+"create_synthetic_reads"+str(rep))
    
    seed_for_wgsim = random.randint(0,1000)

    fastq1 = reference.replace(".fasta", "_1.fastq")
    fastq2 = reference.replace(".fasta", "_2.fastq")
    report = reference.replace(".fasta", "_wgsim_report.txt")

    command = ["wgsim", "-1", "101", "-2", "101", "-d", "300", "-N", str(num_pairs), "-S", str(seed_for_wgsim), "-e", "0.01", "-h", reference, fastq1, fastq2]
    run_command_stdout(command, report)

    return fastq1, fastq2

def run_mcclintock(fastq1, fastq2, reference, consensus, locations, taxonomy, rep, threads, out, config, run_id="", reverse=False):
    if not os.path.exists(out+"/results"):
        os.mkdir(out+"/results")
    
    if not reverse:
        if not os.path.exists(out+"/results/forward"):
            os.mkdir(out+"/results/forward")
        mcc_out = out+"/results/forward/run"+run_id+"_"+str(rep)
    else:
        if not os.path.exists(out+"/results/reverse"):
            os.mkdir(out+"/results/reverse")
        mcc_out = out+"/results/reverse/run"+run_id+"_"+str(rep)

    if not os.path.exists(mcc_out):
        os.mkdir(mcc_out)
    
    mcc_path = config['mcclintock']['path']
    command = ["python3",mcc_path+"/mcclintock.py", "-r", reference, "-c", consensus, "-1", fastq1, "-2", fastq2, "-p", str(threads), "-o", mcc_out, "-g", locations, "-t", taxonomy, "-m", config['mcclintock']['methods'], "--replace_invalid_symbols"]
    if 'augment' in config['mcclintock'].keys() and config['mcclintock']['augment'] is not None:
        command += ["-a", config['mcclintock']['augment']]
    print("running mcclintock... output:", mcc_out)
    run_command_stdout(command, mcc_out+"/run.stdout", log=mcc_out+"/run.stderr")
    if not os.path.exists(mcc_out+"/results/summary/summary_report.txt"):
        sys.stderr.write("run at: "+mcc_out+" failed...")
    

def get_actual_insertions(out):
    actual_insertions = {}
    for f in os.listdir(out):
        if ".bed" in f:
            bed = out+"/"+f
            rep = bed.split("/")[-1]
            rep = rep.split(".")[0]
            rep = "run_"+rep
            with open(bed,"r") as b:
                for line in b:
                    insertion = Insertion()
                    line = line.replace("\n","")
                    split_line = line.split("\t")
                    insertion.chromosome = split_line[0]
                    insertion.start = int(split_line[1])
                    insertion.end = int(split_line[2])
                    insertion.family = split_line[3]
                    # insertion.strand = split_line[5]
            
            actual_insertions[rep] = insertion
    
    return actual_insertions

def get_predicted_insertions(out):
    predicted_insertions = {}
    methods = []
    for d in os.listdir(out):
        rep = d
        rep_num = d.replace("run_","")
        if os.path.exists(out+"/"+d+"/"+rep_num+"/results/"):
            for e in os.listdir(out+"/"+d+"/"+rep_num+"/results/"):
                method = e
                for f in os.listdir(out+"/"+d+"/"+rep_num+"/results/"+e):
                    if "nonredundant.bed" in f:
                        if method not in methods:
                            methods.append(method)
                        with open(out+"/"+d+"/"+rep_num+"/results/"+e+"/"+f, "r") as bed:
                            if method not in predicted_insertions.keys():
                                predicted_insertions[method] = {rep: []}
                            elif rep not in predicted_insertions[method].keys():
                                predicted_insertions[method][rep] = []

                            for line in bed:
                                if "#" not in line:
                                    insertion = Insertion()
                                    line = line.replace("\n","")
                                    split_line = line.split("\t")
                                    if len(split_line) > 5:
                                        insertion.chromosome = split_line[0]
                                        insertion.start = int(split_line[1])
                                        insertion.end = int(split_line[2])
                                        info = split_line[3].split("|")
                                        for x, inf in enumerate(info):
                                            if "reference" in inf:
                                                family = "_".join(info[:x])
                                        insertion.family = family
                                        insertion.strand = split_line[5]
                                        if "non-reference" in line:
                                            insertion.reference = False
                                        else:
                                            insertion.reference = True
                                        
                                        predicted_insertions[method][rep].append(insertion)
    
    return predicted_insertions, methods

def make_out_table(actual_insertions, predicted_insertions, methods, out_csv):
    columns = [["Method","Reference","Nonreference", "Exact", "Within-5", "Within-100", "Within-300", "Within-500"]]
    for method in methods:
        ref_counts = []
        nonref_counts = []
        exacts = []
        within_5s = []
        within_100s = []
        within_300s = []
        within_500s = []
        for rep in predicted_insertions[method].keys():
            ref_count = 0
            nonref_count = 0
            exact = 0
            within_5 = 0
            within_100 = 0
            within_300 = 0
            within_500 = 0
            for insert in predicted_insertions[method][rep]:
                if insert.reference:
                    ref_count +=1
                else:
                    nonref_count +=1
                    if insert.family == actual_insertions[rep].family and insert.chromosome == actual_insertions[rep].chromosome:
                        if (actual_insertions[rep].start == insert.start and insert.end == actual_insertions[rep].end):
                            exact += 1
                        
                        window = 5
                        if ((actual_insertions[rep].start-window <= insert.start and insert.start <= actual_insertions[rep].end+window) or
                        (insert.start <= actual_insertions[rep].start-window and actual_insertions[rep].start-window <= insert.end)):
                            within_5 += 1

                        window = 100
                        if ((actual_insertions[rep].start-window <= insert.start and insert.start <= actual_insertions[rep].end+window) or
                        (insert.start <= actual_insertions[rep].start-window and actual_insertions[rep].start-window <= insert.end)):
                            within_100 += 1

                        window = 300
                        if ((actual_insertions[rep].start-window <= insert.start and insert.start <= actual_insertions[rep].end+window) or
                        (insert.start <= actual_insertions[rep].start-window and actual_insertions[rep].start-window <= insert.end)):
                            within_300 += 1

                        window = 500
                        if ((actual_insertions[rep].start-window <= insert.start and insert.start <= actual_insertions[rep].end+window) or
                        (insert.start <= actual_insertions[rep].start-window and actual_insertions[rep].start-window <= insert.end)):
                            within_500 += 1
            
            ref_counts.append(ref_count)
            nonref_counts.append(nonref_count)
            exacts.append(exact)
            within_5s.append(within_5)
            within_100s.append(within_100)
            within_300s.append(within_300)
            within_500s.append(within_500)
        

        ref_mean = statistics.mean(ref_counts)
        nonref_mean = statistics.mean(nonref_counts)
        exact_mean = statistics.mean(exacts)
        mean_5 = statistics.mean(within_5s)
        mean_100 = statistics.mean(within_100s)
        mean_300 = statistics.mean(within_300s)
        mean_500 = statistics.mean(within_500s)
        columns.append([method, str(ref_mean), str(nonref_mean), str(exact_mean), str(mean_5), str(mean_100), str(mean_300), str(mean_500)])
    
    with open(out_csv, "w") as csv:
        for y in range(0, len(columns[0])):
            line = []
            for x in range(0, len(columns)):
                line.append(columns[x][y])
            line = ",".join(line)
            csv.write(line+"\n")


if __name__ == "__main__":                
    main()