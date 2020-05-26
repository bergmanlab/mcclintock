import argparse
import os
import sys
import traceback
import subprocess
import statistics

# class Value:
#     def __init__(self):
#         self.coverage = ""
#         self.run_name = ""
#         self.run_type = ""
#         self.method = ""
#         self.prediction_type = ""
#         self.prediction = 0

def main():
    args = parse_args()
    if not os.path.exists(args.out+"/summary/"):
        os.mkdir(args.out+"/summary/")

    ref_te_counts = get_ref_te_counts(args.locations, args.taxonomy)

    mcclintock_values, coverages, run_types, methods, pred_types = read_mcclintock_summary(args.out)
    write_allrun_tables(mcclintock_values, coverages, run_types, methods, pred_types, args.out+"/summary/")
    write_summary_tables(mcclintock_values, coverages, run_types, methods, pred_types, args.out+"/summary/")

    


def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock', description="Meta-pipeline to identify transposable element insertions using next generation sequencing data")

    ## required ##
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required=True)
    parser.add_argument("-o", "--out", type=str, help="Directory with output from reference simulation runs [default = '.']", required=True)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required=True)
    

    args = parser.parse_args()

    #check -g
    args.locations = get_abs_path(args.locations)
    #check -t
    args.taxonomy = get_abs_path(args.taxonomy)

    args.out = os.path.abspath(args.out)
    if not os.path.exists(args.out):
        print("-o/--out must point to directory with mcclintock reference simulation results")
        sys.exit(1)
    
    return args


def get_ref_te_counts(tegff, taxontsv):
    te_to_fam = {}
    with open(taxontsv,"r") as tsv:
        for line in tsv:
            line = line.replace("\n","")
            split_line = line.split("\t")
            te_to_fam[split_line[0]] = split_line[1]
    
    fam_counts = {}
    with open(tegff,"r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                te_id = split_line[8].replace("\n","")
                te_id = te_id.split("=")[1]

                fam = te_to_fam[te_id]
                if fam not in fam_counts.keys():
                    fam_counts[fam] = 1
                else:
                    fam_counts[fam] += 1
    
    return fam_counts

def read_mcclintock_summary(out):
    coverages = ["10X", "100X"]
    all_run_types = []
    all_methods = []
    all_pred_types = []
    values = {}

    for cov in coverages:
        reps_dir = out+"/"+cov
        for rep in os.listdir(reps_dir):
            rep_dir = reps_dir+"/"+rep
            for run_type in os.listdir(rep_dir):
                if run_type not in all_run_types:
                    all_run_types.append(run_type)
                csv = rep_dir+"/"+run_type+"/results/summary/te_summary.csv"
                if os.path.exists(csv):
                    with open(csv,"r") as c:
                        for x, line in enumerate(c):
                            line = line.replace("\n","")
                            if x == 0:
                                header = line.split(",")
                                methods = []
                                prediction_types = []
                                for j, r in enumerate(header):
                                    if j > 0:
                                        method = r.split("_")[:-1]
                                        method = "_".join(method)
                                        methods.append(method)
                                        if method not in all_methods:
                                            all_methods.append(method)
                                        prediction_type = r.split("_")[-1]
                                        prediction_types.append(prediction_type)
                                        if prediction_type not in all_pred_types:
                                            all_pred_types.append(prediction_type)
                                totals = [0]*len(methods)
                            else:
                                split_line = line.split(",")
                                for i in range(1, len(split_line)):
                                    val = int(split_line[i])
                                    totals[i-1] += val
                    
                    for col in range(0,len(methods)):
                        if cov not in values.keys():
                            values[cov] = {
                                run_type : { 
                                    prediction_types[col]:{ 
                                        methods[col]: [totals[col]]
                                    } 
                                }
                            }
                        
                        elif run_type not in values[cov].keys():
                            values[cov][run_type] = {
                                prediction_types[col] : {
                                    methods[col] : [totals[col]]
                                }
                            }
                        
                        elif prediction_types[col] not in values[cov][run_type].keys():
                            values[cov][run_type][prediction_types[col]] = {
                                methods[col] : [totals[col]]
                            }
                        
                        elif methods[col] not in values[cov][run_type][prediction_types[col]].keys():
                            values[cov][run_type][prediction_types[col]][methods[col]] = [totals[col]]
                        
                        else:
                            values[cov][run_type][prediction_types[col]][methods[col]].append(totals[col])


    return values, coverages, all_run_types, all_methods, all_pred_types
                                    

def write_allrun_tables(values, coverages, run_types, methods, pred_types, out):
    for pred_type in pred_types:
        out_csv = out+"/"+pred_type+"_runs.csv"        
        for cov in coverages:
            for run_type in run_types:
                out_csv = out+"/"+pred_type+"_"+cov+"_"+run_type+"_runs.csv"
                with open(out_csv,"w") as csv:
                    csv.write(",".join(methods)+"\n")
                    for method in methods:
                        num_runs = len(values[cov][run_type][pred_type][method])
                    for x in range(0,num_runs):
                        out_line = []
                        for method in methods:
                            out_line.append(str(values[cov][run_type][pred_type][method][x]))
                        csv.write(",".join(out_line)+"\n")
                    

def write_summary_tables(values, coverages, run_types, methods, pred_types, out):
    for pred_type in pred_types:
        out_csv = out+"/"+pred_type+"_mean.csv"
        with open(out_csv,"w") as csv:
            csv.write(",".join(["cov","ref_type"]+methods)+"\n")
            for cov in coverages:
                for run_type in run_types:
                    line = [cov, run_type]
                    for method in methods:
                        mean = statistics.mean(values[cov][run_type][pred_type][method])
                        line.append(str(mean))
                    
                    csv.write(",".join(line)+"\n")


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



if __name__ == "__main__":                
    main()