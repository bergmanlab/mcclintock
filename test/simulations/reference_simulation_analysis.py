import argparse
import os
import sys
import traceback
import subprocess
import statistics

def main():
    args = parse_args()
    if not os.path.exists(args.out+"/summary/"):
        os.mkdir(args.out+"/summary/")

    ref_te_counts = get_ref_te_counts(args.locations, args.taxonomy)

    runs = ["default", "ref", "cons", "ref_cons"]
    for run in runs:
        results, methods, families = read_mcclintock_summary( args.out+"/100X", run_type=run)
        write_res_table(results, methods, families, ref_te_counts, args.out+"/summary/100X_"+run+"_results.csv")

    for run in runs:
        results, methods, families = read_mcclintock_summary( args.out+"/10X", run_type=run)
        write_res_table(results, methods, families, ref_te_counts, args.out+"/summary/10X_"+run+"_results.csv")


    


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

def read_mcclintock_summary(out, run_type="default"):
    results = {}
    methods = []
    families = []

    for d in os.listdir(out):
        d_dir = out+"/"+d
        for e in os.listdir(d_dir):
            if e == run_type:
                csv = d_dir+"/"+e+"/results/summary/te_summary.csv"
                if os.path.exists(csv):
                    with open(csv,"r") as c:
                        for x, line in enumerate(c):
                            line = line.replace("\n","")
                            if x == 0:
                                res = line.split(",")
                                for j, r in enumerate(res):
                                    if j > 0 and r not in methods and "_ref" in r:
                                        methods.append(r)
                            
                            else:
                                split_line = line.split(",")
                                if split_line[0] not in families:
                                    families.append(split_line[0])
                                
                                for i in range(1, len(split_line)):
                                    method = res[i]
                                    if "_ref" in method:
                                        family = split_line[0]
                                        val = int(split_line[i])

                                        if method not in results.keys():
                                            results[method] = {family: [val]}
                                        elif family not in results[method].keys():
                                            results[method][family] = [val]
                                        else:
                                            results[method][family].append(val)

    return results, methods, families
                                    


def write_res_table(data, methods, families, ref_tes, out):
    with open(out,"w") as o:
        header = ["TE Family"] + methods
        o.write(",".join(header)+",Actual\n")
        ref_total = 0
        for family in families:
            line = [family]
            for method in methods:
                val = statistics.mean(data[method][family])
                line.append(str(val))

            ref_total += ref_tes[family]
            line.append(str(ref_tes[family]))
            o.write(",".join(line)+"\n")
        
        line = ["total"]
        for method in methods:
            total = 0
            for family in families:
                total += statistics.mean(data[method][family])
            line.append(str(total))
        
        line.append(str(ref_total))
        o.write(",".join(line)+"\n")


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