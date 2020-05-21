import argparse
import os
import sys
import traceback
import subprocess
from multiprocessing import Process, Pool

def main():
    args = parse_args()
    os.mkdir(args.out+"/data")
    fastqs_100x = threaded_make_fastqs(args.reference, args.out, threads=args.proc, start=args.start, end=args.end)
    fastqs_10x = threaded_subsample(args.out, threads=args.proc, start=args.start, end=args.end)
    threaded_mcclintock_run(fastqs_100x, fastqs_10x, args.reference, args.consensus, args.locations, args.taxonomy, args.out, threads=args.proc, start=args.start)




def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock', description="Script to run McClintock on genomes with synthetic reads to measure recovery of reference TEs in each tool")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required=True)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required=True)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required=True)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required=True)
    

    ## optional ##
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required=False)
    parser.add_argument("-s", "--start", type=int, help="the start of the range of seeds to run (default=1)", required=False)
    parser.add_argument("-e", "--end", type=int, help="the end of the range of seeds to run (default=100)", required=False)

    args = parser.parse_args()


    #check -r
    args.reference = get_abs_path(args.reference)
    #check -c
    args.consensus = get_abs_path(args.consensus)
    # check -g
    args.locations = get_abs_path(args.locations)
    # check -t
    args.taxonomy = get_abs_path(args.taxonomy)

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

    if args.start is None:
        args.start = 1
    if args.end is None:
        args.end = 100

    if args.start >= args.end:
        print("-s/--start must be lower than -e/--end")
        sys.exit(1)
    
    return args



def threaded_make_fastqs(ref, out, threads=1, start=1, end=100):
    print("Creating simulated fastq files...")
    fastqs = []
    inputs = []
    for x in range(start, end+1):
        inputs.append([ref, x, out])
        fastqs.append(out+"/data/100X"+str(x)+"simulation")
    
    pool = Pool(processes=threads)
    pool.map(make_fastq, inputs)
    pool.close()
    pool.join()

    return fastqs



def make_fastq(args):
    ref = args[0]
    idx = args[1]
    out = args[2]

    command = ["wgsim", "-1", "101", "-2", "101", "-d", "300", "-N", "6021285", "-S", str(idx), "-e", "0.01", "-h", ref, out+"/data/100X"+str(idx)+"simulation_1.fastq", out+"/data/100X"+str(idx)+"simulation_2.fastq"]
    run_command_stdout(command, out+"/data/wgsim_mutation_report_"+str(idx)+".txt")


def threaded_subsample(out, threads=1, start=1, end=100):
    print("Subsampling fastq files to 10X coverage")
    inputs = []
    fastqs = []
    
    for x in range(start, end+1):
        inputs.append([x, out])
        fastqs.append(out+"/data/10X"+str(x)+"simulation")
    
    pool = Pool(processes=threads)
    pool.map(subsample_fastqs, inputs)
    pool.close()
    pool.join()

    return fastqs


def subsample_fastqs(args):
    idx = args[0]
    out = args[1]
    print("subsampling ", out+"/data/100X"+str(idx)+"simulation_1.fastq")
    command = ["seqtk", "seq", "-s", str(idx), "-f", "0.1", out+"/data/100X"+str(idx)+"simulation_1.fastq"]
    run_command_stdout(command, out+"/data/10X"+str(idx)+"simulation_1.fastq")

    print("subsampling ", out+"/data/100X"+str(idx)+"simulation_2.fastq")
    command = ["seqtk", "seq", "-s", str(idx), "-f", "0.1", out+"/data/100X"+str(idx)+"simulation_2.fastq"]
    run_command_stdout(command, out+"/data/10X"+str(idx)+"simulation_2.fastq")


def threaded_mcclintock_run(fastqs_100x, fastqs_10x, ref, consensus, locations, taxonomy, out, threads=1, start=1):
    print("Starting McClintock runs on simulated reads")
    inputs = []
    os.mkdir(out+"/100X")
    os.mkdir(out+"/10X")
    for x in range(0, len(fastqs_100x)):
        idx = x+start
        os.mkdir(out+"/100X/100X"+str(idx))
        inputs.append([fastqs_100x[x]+"_1.fastq", fastqs_100x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/100X/100X"+str(idx)+"/default", False, False])
        inputs.append([fastqs_100x[x]+"_1.fastq", fastqs_100x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/100X/100X"+str(idx)+"/ref", True, False])
        inputs.append([fastqs_100x[x]+"_1.fastq", fastqs_100x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/100X/100X"+str(idx)+"/cons", False, True])
        inputs.append([fastqs_100x[x]+"_1.fastq", fastqs_100x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/100X/100X"+str(idx)+"/ref_cons", True, True])

        os.mkdir(out+"/10X/10X"+str(idx))
        inputs.append([fastqs_10x[x]+"_1.fastq", fastqs_10x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/10X/10X"+str(idx)+"/default", False, False])
        inputs.append([fastqs_10x[x]+"_1.fastq", fastqs_10x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/10X/10X"+str(idx)+"/ref", True, False])
        inputs.append([fastqs_10x[x]+"_1.fastq", fastqs_10x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/10X/10X"+str(idx)+"/cons", False, True])
        inputs.append([fastqs_10x[x]+"_1.fastq", fastqs_10x[x]+"_2.fastq", ref, consensus, locations, taxonomy, out+"/10X/10X"+str(idx)+"/ref_cons", True, True])
    
    pool = Pool(processes=threads)
    pool.map(mcclintock_run, inputs)
    pool.close()
    pool.join()


def mcclintock_run(args):
    fq1 = args[0]
    fq2 = args[1]
    ref = args[2]
    consensus = args[3]
    locations = args[4]
    taxonomy = args[5]
    out = args[6]
    add_ref = args[7]
    add_cons = args[8]

    os.mkdir(out)

    mcc_path = str(os.path.dirname(os.path.abspath(__file__)))+"/../../"

    command = ["python3",mcc_path+"/mcclintock.py", "-r", ref, "-c", consensus, "-1", fq1, "-2", fq2, "-p", "1", "-o", out, "-g", locations, "-t", taxonomy]

    if add_ref:
        command.append("-R")
    
    if add_cons:
        command.append("-C")

    print("running mcclintock... output:", out)
    run_command_stdout(command, out+"/run.stdout", log=out+"/run.stderr")

    if not os.path.exists(out+"/results/summary/summary_report.txt"):
        sys.stderr.write("run at: "+out+" failed...")






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