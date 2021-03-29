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
    if not os.path.exists(args.out+"/summary/"):
        os.mkdir(args.out+"/summary/")

    actual_insertions = get_actual_insertions(args.out+"/data/forward/")
    predicted_insertions, methods = get_predicted_insertions(args.out+"/results/forward/", args.exclude)
    forward_metrics = make_out_table(actual_insertions, predicted_insertions, methods, args.out+"/summary/forward_summary.csv", args.out+"/summary/forward_metrics.csv")

    actual_insertions = get_actual_insertions(args.out+"/data/reverse/")
    predicted_insertions, methods = get_predicted_insertions(args.out+"/results/reverse/", args.exclude)
    reverse_metrics = make_out_table(actual_insertions, predicted_insertions, methods, args.out+"/summary/reverse_summary.csv", args.out+"/summary/reverse_metrics.csv")

    write_combined_metrics(forward_metrics, reverse_metrics, args.out+"/summary/combined_metrics.csv")


def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock Simulation Analysis', description="Script to summarize the results of a synthetic insertion simulation")

    ## required ##
    parser.add_argument("-o", "--out", type=str, help="Directory containing the simulation results", required=True)

    ## optional ##
    parser.add_argument("-e", "--exclude", type=str, help="bed file of regions in which predictions will be excluded from counts (ex. low recombination regions)", required=False)

    args = parser.parse_args()

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
    
    if args.exclude is not None:
        args.exclude = get_abs_path(args.exclude)

    return args

def get_abs_path(in_file, log=None):
    if os.path.isfile(in_file):
        return os.path.abspath(in_file)
    else:
        msg = " ".join(["Cannot find file:",in_file,"exiting....\n"])
        sys.stderr.write(msg)
        writelog(log, msg)
        sys.exit(1)

def writelog(log, msg):
    if log is not None:
        with open(log, "a") as out:
            out.write(msg)

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

def get_predicted_insertions(out, exclude=None):
    predicted_insertions = {}
    methods = []
    for d in os.listdir(out):
        rep = d
        rep_num = d.replace("run_","")
        if os.path.exists(out+"/"+d+"/"+rep_num+".modref_1/results/"):
            for e in os.listdir(out+"/"+d+"/"+rep_num+".modref_1/results/"):
                method = e
                for f in os.listdir(out+"/"+d+"/"+rep_num+".modref_1/results/"+e):
                    if "nonredundant.bed" in f and "exclude.bed" not in f:
                        if method not in methods:
                            methods.append(method)
                        inbed = out+"/"+d+"/"+rep_num+".modref_1/results/"+e+"/"+f
                        if exclude is not None:
                            tmpbed = out+"/"+d+"/"+rep_num+".modref_1/results/"+e+"/"+f+".exclude.bed"
                            if os.path.exists(tmpbed):
                                subprocess.call(["rm", tmpbed])
                            run_command_stdout(["bedtools", "intersect", "-v", "-a", inbed, "-b", exclude], tmpbed)
                            inbed = tmpbed
                        with open(inbed, "r") as bed:
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

def make_out_table(actual_insertions, predicted_insertions, methods, out_csv, out_metrics):
    columns = [["Method","Reference","Nonreference", "Exact", "Within-5", "Within-100", "Within-300", "Within-500"]]
    rows = [["Method", "TP", "FP","FN", "Precision", "Recall"]]
    for method in methods:
        ref_counts = []
        nonref_counts = []
        exacts = []
        within_5s = []
        within_100s = []
        within_300s = []
        within_500s = []
        all_pred = 0
        true_pos = 0
        false_pos = 0
        false_neg = 0
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
        
            all_pred += nonref_count
            if within_5 > 0:
                true_pos += 1
            else:
                false_neg += 1
            


        ref_mean = statistics.mean(ref_counts)
        nonref_mean = statistics.mean(nonref_counts)
        exact_mean = statistics.mean(exacts)
        mean_5 = statistics.mean(within_5s)
        mean_100 = statistics.mean(within_100s)
        mean_300 = statistics.mean(within_300s)
        mean_500 = statistics.mean(within_500s)
        columns.append([method, str(ref_mean), str(nonref_mean), str(exact_mean), str(mean_5), str(mean_100), str(mean_300), str(mean_500)])

        false_pos = all_pred-true_pos
        if (true_pos+false_pos) > 0:
            precision = true_pos/(true_pos+false_pos)
        else:
            precision = "NaN"
        if (true_pos+false_neg) > 0:
            recall = true_pos/(true_pos+false_neg)
        else:
            recall = "NaN"
        rows.append([method, str(true_pos), str(false_pos), str(false_neg), str(precision), str(recall)])
    
    with open(out_csv, "w") as csv:
        for y in range(0, len(columns[0])):
            line = []
            for x in range(0, len(columns)):
                line.append(columns[x][y])
            line = ",".join(line)
            csv.write(line+"\n")
    
    with open(out_metrics, "w") as csv:
        for row in rows:
            csv.write(",".join(row) + "\n")
    
    return rows

def write_combined_metrics(forward_metrics, reverse_metrics, outfile):
    metrics = {}
    for x,metric in enumerate(forward_metrics):
        if x > 0:
            method = metric[0]
            true_pos = int(metric[1])
            false_pos = int(metric[2])
            false_neg = int(metric[3])
            metrics[method] = [true_pos, false_pos, false_neg]
            
    for x,metric in enumerate(reverse_metrics):
        if x > 0:
            method = metric[0]
            true_pos = int(metric[1])
            false_pos = int(metric[2])
            false_neg = int(metric[3])

            metrics[method] = [metrics[method][0]+true_pos, metrics[method][1]+false_pos, metrics[method][2]+false_neg]
    
    with open(outfile, "w") as out:
        out.write(",".join(forward_metrics[0]) + "\n")
        for method in metrics.keys():
            precision = metrics[method][0] / (metrics[method][0] + metrics[method][1])
            recall = metrics[method][0] / (metrics[method][0] + metrics[method][2])
            outline = [method, str(metrics[method][0]), str(metrics[method][1]), str(metrics[method][2]), str(precision), str(recall)]
            out.write(",".join(outline) + "\n")


if __name__ == "__main__":                
    main()