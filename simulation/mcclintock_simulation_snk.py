#!/usr/bin/env python
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
import glob
import shutil

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

    #check -r
    reference = args.reference
    #check -c
    consensus = args.consensus
    # check -g
    locations = args.locations
    # check -t
    taxonomy = args.taxonomy
    # check -j
    config = args.config

    proc = args.proc
    out = args.out
    start = args.start
    end = args.end
    runid = args.runid
    single = args.single
    coverage = args.coverage
    length = args.length
    insert = args.insert
    error = args.error
    keep_intermediate = args.keep_intermediate
    strand = args.strand
    sim = args.sim
    mcc_version = args.mcc_version
    seed = args.seed

    # # parsing arguments from snakemake
    # #check -r
    # reference = get_abs_path(snakemake.input.ref)
    # #check -c
    # consensus = get_abs_path(snakemake.input.consensus)
    # # check -g
    # locations = get_abs_path(snakemake.input.gff)
    # # check -t
    # taxonomy = get_abs_path(snakemake.input.tax)
    # # check -j
    # config = get_abs_path(snakemake.input.cfg)
    # with open(config, "r") as j:
    #     cfg = json.load(j, object_pairs_hook = OrderedDict)
    # config = cfg

    # #check -p
    # if snakemake.threads is None:
    #     proc = 1
    # else:
    #     proc = int(snakemake.threads)


    # #check -o
    # if snakemake.params.out is None:
    #     out = os.path.abspath(".")
    # else:
    #     out = os.path.abspath(snakemake.params.out)

    #     if not os.path.exists(out):
    #         try:
    #             os.mkdir(out)
    #         except Exception as e:
    #             track = traceback.format_exc()
    #             print(track, file=sys.stderr)
    #             print("cannot create output directory: ",out,"exiting...", file=sys.stderr)
    #             sys.exit(1)

    # # check --start
    # if snakemake.params.start is None:
    #     start = 1
    # else:
    #     start = int(snakemake.params.start)

    # # check --end
    # if snakemake.params.end is None:
    #     end = 1
    # else:
    #     end = int(snakemake.params.end)

    
    # if snakemake.params.runid is None:
    #     runid = ""
    # else:
    #     runid = snakemake.params.runid

    # if snakemake.params.single is None:
    #     single = False
    # else:
    #     single = snakemake.params.single

    # if snakemake.params.coverage is None:
    #     coverage = 100
    # else:
    #     coverage = int(snakemake.params.coverage)

    # if snakemake.params.length is None:
    #     length = 101
    # else:
    #     length = int(snakemake.params.length)

    # if snakemake.params.insertsize is None:
    #     insert = 300
    # else:
    #     insert = int(snakemake.params.insertsize)

    # if snakemake.params.error is None:
    #     error = 0.01
    # else:
    #     error = float(snakemake.params.error)

    # if snakemake.params.keep_intermediate is None:
    #     keep_intermediate = "general"
    # else:
    #     keep_intermediate = snakemake.params.keep_intermediate

    # if snakemake.params.strand is None:
    #     strand = "forward"
    # elif snakemake.params.strand != "forward" and snakemake.params.strand != "reverse":
    #     sys.exit("ERROR: --strand must be forward or reverse \n")
    # else:
    #     strand = snakemake.params.strand


    # if snakemake.params.sim is None:
    #     sim = "wgsim"
    # elif snakemake.params.sim not in ["wgsim", "art"]:
    #     sys.exit("ERROR: --sim must be wgsim or art \n")
    # else:
    #     sim = snakemake.params.sim

    # if snakemake.params.mcc_version is None:
    #     mcc_version = 2
    # else:
    #     mcc_version = int(snakemake.params.mcc_version)
    
    # if snakemake.params.seed is None:
    #     seed = None
    # else:
    #     seed = snakemake.params.seed

    ### start simulation
    for x in range(start,end+1):
        # forward
        consensus_seqs = get_seqs(consensus)
        reference_seqs = get_seqs(reference)

        if strand == "forward":
            reverse = False
            modified_reference = out+"/data/"+str(strand)+"/"+str(runid)+str(x)+".modref.fasta"
            summary_report = out+"/results/"+str(strand)+"/run"+str(runid)+"_"+str(x)+"/"+str(x)+".modref_1/results/summary/te_summary.csv"
            mcc_out = out+"/results/"+str(strand)+"/run"+str(runid)+"_"+str(x)+"/"+str(x)+".modref_1/"
        else:
            reverse = True
            modified_reference = out+"/data/"+str(strand)+"/"+str(runid)+str(x)+".modref.fasta"
            summary_report = out+"/results/"+str(strand)+"/run"+str(runid)+"_"+str(x)+"/"+str(x)+".modref_1/results/summary/te_summary.csv"
            mcc_out = out+"/results/"+str(strand)+"/run"+str(runid)+"_"+str(x)+"/"+str(x)+".modref_1/"

        if not os.path.exists(modified_reference):
            modified_reference = add_synthetic_insertion(reference_seqs, consensus_seqs, config, x, out, run_id=runid, seed=seed, reverse=reverse)

        fastq1 = modified_reference.replace(".fasta", "_1.fastq")
        fastq2 = modified_reference.replace(".fasta", "_2.fastq")
        if not os.path.exists(fastq1) or not os.path.exists(fastq2):
            num_pairs = calculate_num_pairs(modified_reference, coverage, length, single=single)
            fastq1, fastq2 = create_synthetic_reads(sim, modified_reference, coverage, num_pairs, length, insert, error, x, out, run_id=runid, seed=seed)

        if not os.path.exists(summary_report):
            if os.path.exists(mcc_out):
                shutil.rmtree(mcc_out)
            run_mcclintock(fastq1, fastq2, reference, consensus, locations, taxonomy, x, proc, out, config, keep_intermediate, run_id=runid, single=single, reverse=reverse, mcc_version=mcc_version)

def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock Simulation', description="Script to run synthetic insertion simulations to evaluate McClintock component methods")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required=True)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required=True)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required=True)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required=True)
    parser.add_argument("-j","--config", type=str, help="A json config file containing information on TE family TSD size and target sites", required=True)
    ## optional ##
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required=False)
    parser.add_argument("-C","--coverage", type=int, help="The target genome coverage for the simulated reads [default = 100]", required=False)
    parser.add_argument("-l","--length", type=int, help="The read length of the simulated reads [default = 101]", required=False)
    parser.add_argument("-i","--insert", type=int, help="The median insert size of the simulated reads [default = 300]", required=False)
    parser.add_argument("-e","--error", type=float, help="The base error rate for the simulated reads [default = 0.01]", required=False)
    parser.add_argument("-k","--keep_intermediate", type=str, help="This option determines which intermediate files are preserved after McClintock completes [default: general][options: minimal, general, methods, <list,of,methods>, all]", required=False)
    parser.add_argument("--strand", type=str, help="The strand to insert the TE into [options=plus,minus][default = plus]", required=False)
    parser.add_argument("--start", type=int, help="The number of replicates to run. [default = 1]", required=False)
    parser.add_argument("--end", type=int, help="The number of replicates to run. [default = 300]", required=False)
    parser.add_argument("--seed", type=str, help="a seed to the random number generator so runs can be replicated", required=False)
    parser.add_argument("--runid", type=str, help="a string to prepend to output files so that multiple runs can be run at the same time without file name clashes", required=False)
    parser.add_argument("--sim", type=str, help="Short read simulator to use (options=wgsim,art) [default = wgsim]", required=False)
    parser.add_argument("-s","--single", action="store_true", help="runs the simulation in single ended mode", required=False)
    parser.add_argument("--mcc_version", type=int, help="Which version of McClintock to use for the simulation(1 or 2). [default = 2]", required=False, default=2)

    args = parser.parse_args()

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
                os.makedirs(args.out, exist_ok=True)
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

    if args.single is None:
        args.single = False

    if args.coverage is None:
        args.coverage = 100

    if args.length is None:
        args.length = 101

    if args.insert is None:
        args.insert = 300

    if args.error is None:
        args.error = 0.01

    if args.keep_intermediate is None:
        args.keep_intermediate = "general"

    if args.strand is None:
        args.strand = "forward"
    elif args.strand != "forward" and args.strand != "reverse":
        sys.exit("ERROR: --strand must be forward or reverse \n")


    if args.sim is None:
        args.sim = "wgsim"
    elif args.sim not in ["wgsim", "art"]:
        sys.exit("ERROR: --sim must be wgsim or art \n")


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
        os.makedirs(out+"/data", exist_ok=True)

    if not reverse:
        outdir = out+"/data/forward"
    else:
        outdir = out+"/data/reverse"

    if not os.path.exists(outdir):
        os.makedirs(outdir, exist_ok=True)

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


def calculate_num_pairs(fasta, coverage, length, single=False):
    command = ["samtools","faidx", fasta]
    run_command(command)

    total_length = 0
    with open(fasta+".fai", "r") as faidx:
        for line in faidx:
            split_line = line.split("\t")
            contig_length = int(split_line[1])
            total_length += contig_length

    if single:
        num_pairs = (total_length * coverage)/(length)
    else:
        num_pairs = (total_length * coverage)/(2*length)

    return num_pairs

def create_synthetic_reads(simulator, reference, coverage, num_pairs, length, insert, error, rep, out, run_id="", seed=None):
    if seed is not None:
        random.seed(seed+"create_synthetic_reads"+str(rep))
    else:
        random.seed(str(datetime.now())+"create_synthetic_reads"+str(rep))

    seed_for_wgsim = random.randint(0,1000)

    tmp_fastq1 = reference.replace(".fasta", "") + "1.fq"
    tmp_fastq2 = reference.replace(".fasta", "") + "2.fq"

    fastq1 = reference.replace(".fasta", "_1.fastq")
    fastq2 = reference.replace(".fasta", "_2.fastq")
    report = reference.replace(".fasta", "_wgsim_report.txt")

    if simulator == "wgsim":
        command = ["wgsim", "-1", str(length), "-2", str(length), "-d", str(insert), "-N", str(num_pairs), "-S", str(seed_for_wgsim), "-e", str(error), "-h", reference, tmp_fastq1, tmp_fastq2]

    else:
        command = ["art_illumina", "-ss", "HS25", "--rndSeed", str(seed_for_wgsim), "-sam", "-i", reference, "-p", "-l", str(length), "-f", str(coverage), "-m", str(insert), "-s", "10", "-o", reference.replace(".fasta", "")]


    run_command_stdout(command, report)
    run_command(["mv", tmp_fastq1, fastq1])
    run_command(["mv", tmp_fastq2, fastq2])

    return fastq1, fastq2

def run_mcclintock(fastq1, fastq2, reference, consensus, locations, taxonomy, rep, threads, out, config, keep_intermediate, run_id="", reverse=False, single=False, mcc_version=2):
    if not os.path.exists(out+"/results"):
        os.makedirs(out+"/results", exist_ok=True)

    if not reverse:
        if not os.path.exists(out+"/results/forward"):
            os.makedirs(out+"/results/forward", exist_ok=True)
        mcc_out = out+"/results/forward/run"+run_id+"_"+str(rep)
    else:
        if not os.path.exists(out+"/results/reverse"):
            os.makedirs(out+"/results/reverse", exist_ok=True)
        mcc_out = out+"/results/reverse/run"+run_id+"_"+str(rep)

    if not os.path.exists(mcc_out):
        os.makedirs(mcc_out, exist_ok=True)

    if mcc_version == 2:
        mcc_path = config['mcclintock']['path']
        if single:
            command = [
                "python3",mcc_path+"/mcclintock.py",
                    "-r", reference,
                    "-c", consensus,
                    "-1", fastq1,
                    "-p", str(threads),
                    "-o", mcc_out,
                    "-g", locations,
                    "-t", taxonomy,
                    "-m", config['mcclintock']['methods'],
                    "--keep_intermediate", keep_intermediate
            ]
        else:
            command = [
                "python3",mcc_path+"/mcclintock.py",
                    "-r", reference,
                    "-c", consensus,
                    "-1", fastq1,
                    "-2", fastq2,
                    "-p", str(threads),
                    "-o", mcc_out,
                    "-g", locations,
                    "-t", taxonomy,
                    "-m", config['mcclintock']['methods'],
                    "--keep_intermediate", keep_intermediate
            ]

        if 'augment' in config['mcclintock'].keys() and config['mcclintock']['augment'] is not None:
            command += ["-a", config['mcclintock']['augment']]
        print("running mcclintock... output:", mcc_out)
        print(command)
        run_command_stdout(command, mcc_out+"/run.stdout", log=mcc_out+"/run.stderr")
        if not os.path.exists(mcc_out+"/results/summary/te_summary.csv"):
            sys.stderr.write("run at: "+mcc_out+" failed...")
    else:
        mcc_path = config['mcclintock']['v1_path']
        command = [
            mcc_path+"/mcclintock.sh",
            "-o", mcc_out,
            "-r", reference,
            "-c", consensus,
            "-g", locations,
            "-t", taxonomy,
            "-1", fastq1,
            "-p", str(threads),
            "-i"
        ]
        if not single:
            command += ["-2", fastq2]

        if 'augment' in config['mcclintock'].keys() and config['mcclintock']['augment'] is not None:
            command += ["-C"]

        print("running mcclintock... output:", mcc_out)
        print(command)
        run_command_stdout(command, mcc_out+"/run.stdout", log=mcc_out+"/run.stderr")
        reorder_mcc1_output(rep, mcc_out)


def reorder_mcc1_output(rep, out):
    os.makedirs(f"{out}/{rep}.modref_1/", exist_ok=True)
    results_dir = f"{out}/{rep}.modref_1/results/"
    os.makedirs(results_dir, exist_ok=True)
    beds = glob.glob(out+'/*/*/results/*_nonredundant.bed')
    for bed in beds:
        base_name = bed.split("/")[-1]
        method = base_name.replace("_nonredundant.bed","")
        method = method[method.find("_")+1:]
        method_dir = results_dir+"/"+method
        os.makedirs(method_dir, exist_ok=True)
        with open(bed,"r") as inbed, open(method_dir+"/"+base_name, "w") as outbed:
            for line in inbed:
                outbed.write(line.replace("_","|", 1))








if __name__ == "__main__":
    main()
