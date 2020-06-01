import argparse
import os
import sys
import traceback
import subprocess
import statistics


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
    
    actual_insertions = get_actual_insertions(args.out+"/forward_data")
    predicted_insertions, methods = get_predicted_insertions(args.out+"/forward_results")
    make_out_table(actual_insertions, predicted_insertions, methods, args.out+"/summary/forward_summary.csv")

    actual_insertions = get_actual_insertions(args.out+"/reverse_data")
    predicted_insertions, methods = get_predicted_insertions(args.out+"/reverse_results")
    make_out_table(actual_insertions, predicted_insertions, methods, args.out+"/summary/reverse_summary.csv")



def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock', description="Meta-pipeline to identify transposable element insertions using next generation sequencing data")

    ## required ##
    parser.add_argument("-o", "--out", type=str, help="Directory with output from reference simulation runs [default = '.']", required=True)

    args = parser.parse_args()

    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        print("-o/--out must point to directory with mcclintock reference simulation results")
        sys.exit(1)
    
    return args


def get_actual_insertions(out):
    actual_insertions = {}
    for f in os.listdir(out):
        if ".bed" in f:
            bed = out+"/"+f
            name = f.replace(".bed","")
            with open(bed,"r") as b:
                for line in b:
                    insertion = Insertion()
                    line = line.replace("\n","")
                    split_line = line.split("\t")
                    insertion.chromosome = split_line[0]
                    insertion.start = int(split_line[1])
                    insertion.end = int(split_line[2])
                    insertion.family = split_line[3]
                    insertion.strand = split_line[5]
            
            actual_insertions[name] = insertion
    
    return actual_insertions
    

def get_predicted_insertions(out):
    predicted_insertions = {}
    methods = []
    for d in os.listdir(out):
        rep = d
        for e in os.listdir(out+"/"+d+"/results/"):
            method = e
            for f in os.listdir(out+"/"+d+"/results/"+e):
                if "nonredundant.bed" in f:
                    if method not in methods:
                        methods.append(method)
                    with open(out+"/"+d+"/results/"+e+"/"+f, "r") as bed:
                        for line in bed:
                            if "#" not in line:
                                insertion = Insertion()
                                line = line.replace("\n","")
                                split_line = line.split("\t")
                                if len(split_line) > 5:
                                    insertion.chromosome = split_line[0]
                                    insertion.start = int(split_line[1])
                                    insertion.end = int(split_line[2])
                                    info = split_line[3].split("_")
                                    for x, inf in enumerate(info):
                                        if "reference" in inf:
                                            family = "_".join(info[:x])
                                    insertion.family = family
                                    insertion.strand = split_line[5]
                                    if "non-reference" in line:
                                        insertion.reference = False
                                    else:
                                        insertion.reference = True
                                    
                                    if method not in predicted_insertions.keys():
                                        predicted_insertions[method] = {rep: [insertion]}
                                    elif rep not in predicted_insertions[method].keys():
                                        predicted_insertions[method][rep] = [insertion]
                                    else:
                                        predicted_insertions[method][rep].append(insertion)
    
    return predicted_insertions, methods


def make_out_table(actual_insertions, predicted_insertions, methods, out_csv):
    columns = [["Method","Reference","Nonreference", "Exact", "Within-100", "Within-300", "Within-500"]]
    for method in methods:
        ref_counts = []
        nonref_counts = []
        exacts = []
        within_100s = []
        within_300s = []
        within_500s = []
        for rep in predicted_insertions[method].keys():
            ref_count = 0
            nonref_count = 0
            exact = 0
            within_100 = 0
            within_300 = 0
            within_500 = 0
            for insert in predicted_insertions[method][rep]:
                if insert.reference:
                    ref_count +=1
                else:
                    nonref_count +=1
                    if ((insert.family == actual_insertions[rep].family) and (actual_insertions[rep].start == insert.start and insert.start == actual_insertions[rep].end)):
                        exact += 1
                    
                    window = 100
                    if ((insert.family == actual_insertions[rep].family) and ((actual_insertions[rep].start-window <= insert.start and insert.start <= actual_insertions[rep].end+window) or
                       (insert.start <= actual_insertions[rep].start-window and actual_insertions[rep].start-window <= insert.end))):
                        within_100 += 1

                    window = 300
                    if ((insert.family == actual_insertions[rep].family) and ((actual_insertions[rep].start-window <= insert.start and insert.start <= actual_insertions[rep].end+window) or
                       (insert.start <= actual_insertions[rep].start-window and actual_insertions[rep].start-window <= insert.end))):
                        within_300 += 1

                    window = 500
                    if ((insert.family == actual_insertions[rep].family) and ((actual_insertions[rep].start-window <= insert.start and insert.start <= actual_insertions[rep].end+window) or
                       (insert.start <= actual_insertions[rep].start-window and actual_insertions[rep].start-window <= insert.end))):
                        within_500 += 1
            
            ref_counts.append(ref_count)
            nonref_counts.append(nonref_count)
            exacts.append(exact)
            within_100s.append(within_100)
            within_300s.append(within_300)
            within_500s.append(within_500)
        
        ref_mean = statistics.mean(ref_counts)
        nonref_mean = statistics.mean(nonref_counts)
        exact_mean = statistics.mean(exacts)
        mean_100 = statistics.mean(within_100s)
        mean_300 = statistics.mean(within_300s)
        mean_500 = statistics.mean(within_500s)
        columns.append([method, str(ref_mean), str(nonref_mean), str(exact_mean), str(mean_100), str(mean_300), str(mean_500)])
    
    with open(out_csv, "w") as csv:
        for y in range(0, len(columns[0])):
            line = []
            for x in range(0, len(columns)):
                line.append(columns[x][y])
            line = ",".join(line)
            csv.write(line+"\n")




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

if __name__ == "__main__":                
    main()