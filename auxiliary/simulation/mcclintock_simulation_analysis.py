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
    rows = [["Method", "TP_0", "FP_0","FN_0", "Precision_0", "Recall_0", "TP_5", "FP_5","FN_5", "Precision_5", "Recall_5", "TP_100", "FP_100","FN_100", "Precision_100", "Recall_100", "TP_300", "FP_300","FN_300", "Precision_300", "Recall_300", "TP_500", "FP_500","FN_500", "Precision_500", "Recall_500"]]
    for method in methods:
        ref_counts = []
        nonref_counts = []
        exacts = []
        within_5s = []
        within_100s = []
        within_300s = []
        within_500s = []
        all_pred = 0
        true_pos_0 = true_pos_5 = true_pos_100 = true_pos_300 = true_pos_500 = 0
        false_pos_0 = false_pos_5 = false_pos_100 = false_pos_300 = false_pos_500 = 0
        false_neg_0 = false_neg_5 = false_neg_100 = false_neg_300 = false_neg_500 = 0
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
            if exact > 0:
                true_pos_0 += 1
            else:
                false_neg_0 += 1
            if within_5 > 0:
                true_pos_5 += 1
            else:
                false_neg_5 += 1
            if within_100 > 0:
                true_pos_100 += 1
            else:
                false_neg_100 += 1
            if within_300 > 0:
                true_pos_300 += 1
            else:
                false_neg_300 += 1
            if within_500 > 0:
                true_pos_500 += 1
            else:
                false_neg_500 += 1
            


        ref_mean = statistics.mean(ref_counts)
        nonref_mean = statistics.mean(nonref_counts)
        exact_mean = statistics.mean(exacts)
        mean_5 = statistics.mean(within_5s)
        mean_100 = statistics.mean(within_100s)
        mean_300 = statistics.mean(within_300s)
        mean_500 = statistics.mean(within_500s)
        columns.append([method, str(ref_mean), str(nonref_mean), str(exact_mean), str(mean_5), str(mean_100), str(mean_300), str(mean_500)])

        false_pos_0 = all_pred-true_pos_0
        if (true_pos_0+false_pos_0) > 0:
            precision_0 = true_pos_0/(true_pos_0+false_pos_0)
        else:
            precision_0 = "NaN"
        if (true_pos_0+false_neg_0) > 0:
            recall_0 = true_pos_0/(true_pos_0+false_neg_0)
        else:
            recall_0 = "NaN"
        false_pos_5 = all_pred-true_pos_5
        if (true_pos_5+false_pos_5) > 0:
            precision_5 = true_pos_5/(true_pos_5+false_pos_5)
        else:
            precision_5 = "NaN"
        if (true_pos_5+false_neg_5) > 0:
            recall_5 = true_pos_5/(true_pos_5+false_neg_5)
        else:
            recall_5 = "NaN"
        false_pos_100 = all_pred-true_pos_100
        if (true_pos_100+false_pos_100) > 0:
            precision_100 = true_pos_100/(true_pos_100+false_pos_100)
        else:
            precision_100 = "NaN"
        if (true_pos_100+false_neg_100) > 0:
            recall_100 = true_pos_100/(true_pos_100+false_neg_100)
        else:
            recall_100 = "NaN"
        false_pos_300 = all_pred-true_pos_300
        if (true_pos_300+false_pos_300) > 0:
            precision_300 = true_pos_300/(true_pos_300+false_pos_300)
        else:
            precision_300 = "NaN"
        if (true_pos_300+false_neg_300) > 0:
            recall_300 = true_pos_300/(true_pos_300+false_neg_300)
        else:
            recall_300 = "NaN"
        false_pos_500 = all_pred-true_pos_500
        if (true_pos_500+false_pos_500) > 0:
            precision_500 = true_pos_500/(true_pos_500+false_pos_500)
        else:
            precision_500 = "NaN"
        if (true_pos_500+false_neg_500) > 0:
            recall_500 = true_pos_500/(true_pos_500+false_neg_500)
        else:
            recall_500 = "NaN"
        rows.append([method, str(true_pos_0), str(false_pos_0), str(false_neg_0), str(precision_0), str(recall_0), str(true_pos_5), str(false_pos_5), str(false_neg_5), str(precision_5), str(recall_5), str(true_pos_100), str(false_pos_100), str(false_neg_100), str(precision_100), str(recall_100), str(true_pos_300), str(false_pos_300), str(false_neg_300), str(precision_300), str(recall_300), str(true_pos_500), str(false_pos_500), str(false_neg_500), str(precision_500), str(recall_500)])
    
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
            true_pos_0 = int(metric[1])
            false_pos_0 = int(metric[2])
            false_neg_0 = int(metric[3])
            true_pos_5 = int(metric[6])
            false_pos_5 = int(metric[7])
            false_neg_5 = int(metric[8])
            true_pos_100 = int(metric[11])
            false_pos_100 = int(metric[12])
            false_neg_100 = int(metric[13])
            true_pos_300 = int(metric[16])
            false_pos_300 = int(metric[17])
            false_neg_300 = int(metric[18])
            true_pos_500 = int(metric[21])
            false_pos_500 = int(metric[22])
            false_neg_500 = int(metric[23])
            metrics[method] = [true_pos_0, false_pos_0, false_neg_0, true_pos_5, false_pos_5, false_neg_5, true_pos_100, false_pos_100, false_neg_100, true_pos_300, false_pos_300, false_neg_300, true_pos_500, false_pos_500, false_neg_500]
            
    for x,metric in enumerate(reverse_metrics):
        if x > 0:
            method = metric[0]
            true_pos_0 = int(metric[1])
            false_pos_0 = int(metric[2])
            false_neg_0 = int(metric[3])
            true_pos_5 = int(metric[6])
            false_pos_5 = int(metric[7])
            false_neg_5 = int(metric[8])
            true_pos_100 = int(metric[11])
            false_pos_100 = int(metric[12])
            false_neg_100 = int(metric[13])
            true_pos_300 = int(metric[16])
            false_pos_300 = int(metric[17])
            false_neg_300 = int(metric[18])
            true_pos_500 = int(metric[21])
            false_pos_500 = int(metric[22])
            false_neg_500 = int(metric[23])
            metrics[method] = [metrics[method][0]+true_pos_0, metrics[method][1]+false_pos_0, metrics[method][2]+false_neg_0, metrics[method][3]+true_pos_5, metrics[method][4]+false_pos_5, metrics[method][5]+false_neg_5, metrics[method][6]+true_pos_100, metrics[method][7]+false_pos_100, metrics[method][8]+false_neg_100, metrics[method][9]+true_pos_300, metrics[method][10]+false_pos_300, metrics[method][11]+false_neg_300, metrics[method][12]+true_pos_500, metrics[method][13]+false_pos_500, metrics[method][14]+false_neg_500]
    
    with open(outfile, "w") as out:
        out.write(",".join(forward_metrics[0]) + "\n")
        for method in metrics.keys():
            if metrics[method][0] + metrics[method][1] > 0:
                precision_0 = metrics[method][0] / (metrics[method][0] + metrics[method][1])
            else:
                precision_0 = 0
            if metrics[method][0] + metrics[method][2] > 0:
                recall_0 = metrics[method][0] / (metrics[method][0] + metrics[method][2])
            else:
                recall_0 = 0

            if metrics[method][3] + metrics[method][4] > 0:
                precision_5 = metrics[method][3] / (metrics[method][3] + metrics[method][4])
            else:
                precision_5 = 0
            if metrics[method][3] + metrics[method][5] > 0:
                recall_5 = metrics[method][3] / (metrics[method][3] + metrics[method][5])
            else:
                recall_5 = 0

            if metrics[method][6] + metrics[method][7] > 0:
                precision_100 = metrics[method][6] / (metrics[method][6] + metrics[method][7])
            else:
                precision_100 = 0
            if metrics[method][6] + metrics[method][8] > 0:
                recall_100 = metrics[method][6] / (metrics[method][6] + metrics[method][8])
            else:
                recall_100 = 0

            if metrics[method][9] + metrics[method][10] > 0:
                precision_300 = metrics[method][9] / (metrics[method][9] + metrics[method][10])
            else:
                precision_300 = 0
            if metrics[method][9] + metrics[method][11] > 0:
                recall_300 = metrics[method][9] / (metrics[method][9] + metrics[method][11])
            else:
                recall_300 = 0

            if metrics[method][12] + metrics[method][13] > 0:
                precision_500 = metrics[method][12] / (metrics[method][12] + metrics[method][13])
            else:
                precision_500 = 0
            if metrics[method][12] + metrics[method][14] > 0:
                recall_500 = metrics[method][12] / (metrics[method][12] + metrics[method][14])
            else:
                recall_500 = 0

            outline = [method, str(metrics[method][0]), str(metrics[method][1]), str(metrics[method][2]), str(precision_0), str(recall_0), str(metrics[method][3]), str(metrics[method][4]), str(metrics[method][5]), str(precision_5), str(recall_5), str(metrics[method][6]), str(metrics[method][7]), str(metrics[method][8]), str(precision_100), str(recall_100), str(metrics[method][9]), str(metrics[method][10]), str(metrics[method][11]), str(precision_300), str(recall_300), str(metrics[method][12]), str(metrics[method][13]), str(metrics[method][14]), str(precision_500), str(recall_500)]
            out.write(",".join(outline) + "\n")


if __name__ == "__main__":                
    main()