import argparse
import os
import sys
import traceback
import subprocess
from multiprocessing import Process, Pool
from Bio import SeqIO
from Bio.Seq import Seq


class ModifiedChrom:
    def __init__(self):
        self.name = ""
        self.seq = ""
        self.insert = ""
        self.insert_start = -1
        self.insert_end = -1
        self.strand = "?"

class Feature:
    def __init__(self):
        self.chrom = ""
        self.start = -1
        self.end = -1
        self.strand = "?"

def main():
    args = parse_args()
    
    trnas = get_trnas(args.trna)
    chroms_with_inserts = make_inserts(args.reference, args.consensus, trnas, strand="+")
    if not os.path.exists(args.out+"/forward_data"):
        os.mkdir(args.out+"/forward_data")
    fastas = make_fastas(args.reference, chroms_with_inserts, args.out+"/forward_data", start=args.start, end=args.end)
    make_beds(chroms_with_inserts, args.out+"/forward_data", start=args.start, end=args.end)
    fastqs = threaded_make_fastqs(fastas, args.out+"/forward_data", threads=args.proc)
    if not os.path.exists(args.out+"/forward_results"):
        os.mkdir(args.out+"/forward_results")
    threaded_mcclintock_run(fastqs, args.reference, args.consensus, args.locations, args.taxonomy, args.out+"/forward_results", threads=args.proc)



    if not os.path.exists(args.out+"/reverse_data"):
        os.mkdir(args.out+"/reverse_data")
    chroms_with_inserts = make_inserts(args.reference, args.consensus, trnas, strand="-")
    fastas = make_fastas(args.reference, chroms_with_inserts, args.out+"/reverse_data", start=args.start, end=args.end)
    make_beds(chroms_with_inserts, args.out+"/reverse_data", start=args.start, end=args.end)
    fastqs = threaded_make_fastqs(fastas, args.out+"/reverse_data", threads=args.proc)
    if not os.path.exists(args.out+"/reverse_results"):
        os.mkdir(args.out+"/reverse_results")
    threaded_mcclintock_run(fastqs, args.reference, args.consensus, args.locations, args.taxonomy, args.out+"/reverse_results", threads=args.proc)



def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock', description="Meta-pipeline to identify transposable element insertions using next generation sequencing data")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required=True)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required=True)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required=True)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required=True)
    parser.add_argument("--trna", type=str, help="A gff indicating the locations of tRNAs in the reference genome", required=True)
    

    ## optional ##
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required=False)
    parser.add_argument("-s", "--start", type=int, help="the start of the range of seeds to run (default=1)", required=False)
    parser.add_argument("-e", "--end", type=int, help="the end of the range of seeds to run (default=299)", required=False)

    args = parser.parse_args()


    #check -r
    args.reference = get_abs_path(args.reference)
    #check -c
    args.consensus = get_abs_path(args.consensus)
    # check -g
    args.locations = get_abs_path(args.locations)
    # check -t
    args.taxonomy = get_abs_path(args.taxonomy)
    # check --trna
    args.trna = get_abs_path(args.trna)

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
        args.end = 299

    if args.start >= args.end:
        print("-s/--start must be lower than -e/--end")
        sys.exit(1)
    
    return args


def get_trnas(trna_gff):
    trnas = []
    with open(trna_gff, "r") as gff:
        for line in gff:
            if "#" not in line:
                split_line = line.split("\t")
                trna = Feature()
                trna.chrom = split_line[0]
                trna.start = int(split_line[3])
                trna.end = int(split_line[4])
                trna.strand = split_line[6]

                trnas.append(trna)
    
    return trnas

def make_inserts(reference, consensus, trnas, strand="+"):
    ref_records = SeqIO.parse(reference, "fasta")
    ref_chroms = {}
    for record in ref_records:
        ref_chroms[str(record.id)] = str(record.seq)

    consensus_records = SeqIO.parse(consensus, "fasta")
    
    consensus_contigs = {}
    for record in consensus_records:
        seq = Seq(str(record.seq))
        if strand != "+":
            seq = seq.reverse_complement()
        consensus_contigs[str(record.id)] = str(seq)
    
    
    te_types = ["TY1", "TY2", "TY3", "TY4"]
    te_idx = 0
    synthetic_insertions = []
    for trna in trnas:
        if te_idx >= len(te_types):
            te_idx = 0

        insertion = ModifiedChrom()
        seq = ref_chroms[trna.chrom]
        if trna.strand == "+":
            if te_types[te_idx] == "TY3":
                seq_start = seq[0:trna.start-12]
                seq_end = seq[trna.start-17:len(seq)]
            else:
                seq_start = seq[0:trna.start-195]
                seq_end = seq[trna.start-200:len(seq)]
        
        else:
            if te_types[te_idx] == "TY3":
                seq_start = seq[0:trna.start+17]
                seq_end = seq[trna.start+12:len(seq)]
            else:
                seq_start = seq[0:trna.start+200]
                seq_end = seq[trna.start+195:len(seq)]
                # print(seq_start[-5:], seq_end[0:5], len(seq_start+seq_end) - len(seq))
        
        
        insertion.strand = strand
        insertion.name = trna.chrom
        insertion.insert = te_types[te_idx]
        insertion.seq = seq_start+consensus_contigs[te_types[te_idx]]+seq_end
        insertion.insert_start = len(seq_start)
        insertion.insert_end = len(seq_start+consensus_contigs[te_types[te_idx]])
        synthetic_insertions.append(insertion)
        te_idx += 1

    return synthetic_insertions


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

def make_fastas(reference, chroms_with_inserts, out, start=1, end=299):
    fastas = []
    for x in range(start-1, end):
        out_fasta = out+"/insertion"+str(x+1)+"_genome.fasta"
        with open(out_fasta, "w") as outfa:
            ref_records = SeqIO.parse(reference, "fasta")
            for record in ref_records:
                if str(record.id) == chroms_with_inserts[x].name:
                    print("replacing chrom with chrom with synthetic insert", x+1)
                    seq = chroms_with_inserts[x].seq
                else:
                    seq = str(record.seq)
                
                outfa.write(">"+str(record.id)+"\n")
                outfa.write(seq+"\n")
        
        fastas.append(out_fasta)
    
    return fastas

def make_beds(chroms_with_inserts, out, start=1, end=299):
    for x in range(start-1, end):
        out_bed = out+"/insertion"+str(x+1)+"_genome.bed"
        with open(out_bed,"w") as bed:
            line = [chroms_with_inserts[x].name, str(chroms_with_inserts[x].insert_start), str(chroms_with_inserts[x].insert_end), chroms_with_inserts[x].insert, "0", chroms_with_inserts[x].strand]
            bed.write("\t".join(line)+"\n")

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


def threaded_make_fastqs(refs, out, threads=1):
    print("Creating simulated fastq files...")
    fastqs = []
    inputs = []
    for ref in refs:
        num_pairs = calculate_num_pairs(ref)
        fastq_base_name = ref.replace(".fasta", "")
        inputs.append([ref, num_pairs, fastq_base_name])
        fastqs.append(fastq_base_name)
    
    pool = Pool(processes=threads)
    pool.map(make_fastq, inputs)
    pool.close()
    pool.join()

    return fastqs



def make_fastq(args):
    ref = args[0]
    num_pairs = args[1]
    fq_base_name = args[2]

    command = ["wgsim", "-1", "101", "-2", "101", "-d", "300", "-N", str(num_pairs), "-S", "42", "-e", "0.01", "-h", ref, fq_base_name+"_1.fastq", fq_base_name+"_2.fastq"]
    run_command_stdout(command, fq_base_name+"_wgsim_report.txt")


def threaded_mcclintock_run(fastqs, ref, consensus, locations, taxonomy, out, threads=1):
    print("Starting McClintock runs on simulated reads")
    inputs = []
    for fq in fastqs:
        basename = fq.split("/")[-1]
        os.mkdir(out+"/"+basename)
        inputs.append([fq+"_1.fastq", fq+"_2.fastq", ref, consensus, locations, taxonomy, out+"/"+basename, False, False])
    
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

if __name__ == "__main__":                
    main()