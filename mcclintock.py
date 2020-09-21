#!/usr/bin/env python3

import argparse
import os
import sys
import json
import random
import gzip
from datetime import datetime
import traceback

try:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    import scripts.mccutils as mccutils
    import config.config as config
    import config.install.install as config_install
    from Bio import SeqIO
except ImportError as e:
    print(e)
    sys.exit("ERROR...unable to load required python modules\ntry reinstalling and/or activating mcclintock environment:\n\tconda env create -f install/envs/mcclintock.yml --name mcclintock\n\tconda activate mcclintock\n")





def main():
    full_command = " ".join(["python3"] + sys.argv)
    current_directory = os.getcwd()
    args = parse_args()
    check_input_files(args.reference, args.consensus, args.first, fq2=args.second, locations=args.locations, taxonomy=args.taxonomy, coverage_fasta=args.coverage_fasta, augment_fasta=args.augment)
    mccutils.mkdir(args.out+"/logs")
    mccutils.mkdir(args.out+"/tmp")
    sample_name = mccutils.get_base_name(args.first, fastq=True)
    ref_name = mccutils.get_base_name(args.reference)
    run_id = make_run_config(args, sample_name, ref_name, full_command, current_directory)
    run_workflow(args, sample_name, ref_name, run_id, debug=args.debug, annotations_only=args.make_annotations)

def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock', description="Meta-pipeline to identify transposable element insertions using next generation sequencing data")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required='--install' not in sys.argv)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required='--install' not in sys.argv)
    parser.add_argument("-1", "--first", type=str, help="The path of the first fastq file from paired end read sequencing or the fastq file from single read sequencing", required='--install' not in sys.argv)
    

    ## optional ##
    parser.add_argument("-2", "--second", type=str, help="The path of the second fastq file from a paired end read sequencing", required=False)
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required=False)
    parser.add_argument("-m", "--methods", type=str, help="A comma-delimited list containing the software you want the pipeline to use for analysis. e.g. '-m relocate,TEMP,ngs_te_mapper' will launch only those three methods", required=False)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required=False)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required=False)
    parser.add_argument("-s", "--coverage_fasta", type=str, help="A fasta file that will be used for TE-based coverage analysis, if not supplied then the consensus sequences of the TEs will be used for the analysis", required=False)
    parser.add_argument("-T", "--comments", action="store_true", help="If this option is specified then fastq comments (e.g. barcode) will be incorporated to SAM output. Warning: do not use this option if the input fastq files do not have comments", required=False)
    # parser.add_argument("-b", "--keep_bam", action="store_true", help="Retain the sorted and indexed BAM file of the paired end data aligned to the reference genome", required=False)
    # parser.add_argument("-i", "--remove_intermediate", action="store_true", help="If this option is specified then all sample specific intermediate files will be removed, leaving only the overall results. The default is to leave sample specific intermediate files", required=False)
    parser.add_argument("-a", "--augment", type=str, help="A fasta file of TE sequences that will be included as extra chromosomes in the reference file (useful if the organism is known to have TEs that are not present in the reference strain)", required=False)
    parser.add_argument("--clean", action="store_true", help="This option will make sure mcclintock runs from scratch and doesn't reuse files already created", required=False)
    parser.add_argument("--install", action="store_true", help="This option will install the dependencies of mcclintock", required=False)
    parser.add_argument("--debug", action="store_true", help="This option will allow snakemake to print progress to stdout", required=False)
    parser.add_argument("--slow", action="store_true", help="This option runs without attempting to optimize thread usage to run rules concurrently. Each multithread rule will use the max processors designated by -p/--proc", required=False)
    parser.add_argument("--make_annotations", action="store_true", help="This option will only run the pipeline up to the creation of the repeat annotations", required=False)

    args = parser.parse_args()

    if args.debug is None:
        args.debug = False

    if args.install:
        mccutils.log("install","installing dependencies")
        mccutils.log("install","WARNING: this could take awhile")
        install(clean=args.clean, debug=args.debug)
        sys.exit(0)

    #check -r
    args.reference = mccutils.get_abs_path(args.reference)
    #check -c
    args.consensus = mccutils.get_abs_path(args.consensus)
    #check -1
    args.first = mccutils.get_abs_path(args.first)
    #check -2
    if args.second is not None:
        args.second = mccutils.get_abs_path(args.second)

    #check -p
    if args.proc is None:
        args.proc = 1

    #check -o
    if args.out is None:
        args.out = os.path.abspath(".")
    else:
        args.out = os.path.abspath(args.out)
        try:
            mccutils.mkdir(args.out)
        except Exception as e:
            track = traceback.format_exc()
            print(track, file=sys.stderr)
            print("cannot create output directory: ",args.out,"exiting...", file=sys.stderr)
            sys.exit(1)

    

    #check -m
    # If only one fastq has been supplied assume this is single ended data and launch only ngs_te_mapper and RelocaTE
    if args.second is None:
        valid_methods = config.SINGLE_END_METHODS #from config.py
    else:
        valid_methods = config.ALL_METHODS #from config.py
    
    if args.methods is None:
        args.methods = valid_methods
    
    else:
        args.methods = args.methods.split(",")
        for x,method in enumerate(args.methods):
            args.methods[x] = method.lower()
            if args.methods[x] not in valid_methods:
                sys.stderr.write(" ".join(["Method:",method, "not a valid method...", "Valid methods:"," ".join(valid_methods),"\n"]))
                sys.exit(1)

    # check -g
    if args.locations is not None:
        args.locations = mccutils.get_abs_path(args.locations)

        if args.taxonomy is None:
            sys.stderr.write("If a GFF file is supplied (-g/--locations) then a TE taxonomy file that links it to the fasta consensus is also needed (-t/--taxonomy)...exiting...\n")
            sys.exit(1)
    
    # check -t
    if args.taxonomy is not None:
        args.taxonomy = mccutils.get_abs_path(args.taxonomy)
    

    # check -s
    if args.coverage_fasta is not None:
        args.coverage_fasta = mccutils.get_abs_path(args.coverage_fasta)

    # check -T
    if args.comments is None:
        args.comments = False

    # check -a
    if args.augment is not None:
        args.augment = mccutils.get_abs_path(args.augment)
    
    return args

def check_input_files(ref, consensus, fq1, fq2=None, locations=None, taxonomy=None, coverage_fasta=None, augment_fasta=None):
    # check reference fasta
    mccutils.log("setup","checking reference fasta: "+ref)
    try:
        with open(ref,"r") as fa:
            for record in SeqIO.parse(fa, "fasta"):
                pass
    except Exception as e:
        print(e)
        sys.exit(ref+" appears to be a malformed FastA file..exiting...\n")

    mccutils.log("setup","checking fq1: "+fq1)
    #check fq1
    if ".fastq" not in fq1 and ".fq" not in fq1:
        sys.exit(fq1+" is not a (.fastq/.fq) file, exiting...\n")
    
    #check fq2
    if fq2 is not None:
        mccutils.log("setup","checking fq2: "+fq2)
        if ".fastq" not in fq2 and ".fq" not in fq2:
            sys.exit(fq2+" is not a (.fastq/.fq) file, exiting...\n")
    

    # checking consensus
    mccutils.log("setup","checking consensus fasta: "+consensus)
    consensus_seq_names = []
    try:
        with open(consensus,"r") as fa:
            for record in SeqIO.parse(fa, "fasta"):
                consensus_seq_names.append(str(record.id))
    except Exception as e:
        print(e)
        sys.exit(consensus+" appears to be a malformed FastA file..exiting...\n")
    
    #check locations gff
    gff_ids = []
    if locations is not None:
        mccutils.log("setup","checking locations gff: "+locations)
        with open(locations,"r") as gff:
            for line in gff:
                if "#" not in line:
                    split_line = line.split("\t")
                    if len(split_line) < 9:
                        sys.exit(locations+" appears to be a malformed GFF file..exiting...\n")
                    else:
                        feats = split_line[8]
                        split_feats = feats.split(";")
                        gff_id = ""
                        for feat in split_feats:
                            if feat[:3] == "ID=":
                                gff_id = feat.split("=")[1].replace("\n","")
                                if gff_id not in gff_ids:
                                    gff_ids.append(gff_id)
                                else:
                                    sys.exit("ID: "+gff_id+" is not unique. please ensure each feature has a unique ID\n")
                        if gff_id == "":
                            sys.exit("GFF line: "+line+" is missing an ID attribute (ex. ID=chr1_TY1s1+\n")
                        
    
    # check taxonomy
    if taxonomy is not None:
        mccutils.log("setup","checking taxonomy TSV: "+taxonomy)
        with open(taxonomy, "r") as tsv:
            for line in tsv:
                split_line = line.split("\t")
                if len(split_line) != 2:
                    sys.exit(taxonomy+" does not have two columns. Should be tab-separated file with feature ID and TE family as columns\n")
                else:
                    te_id = split_line[0]
                    te_family = split_line[1].replace("\n","")

                    if te_id not in gff_ids:
                        sys.exit("TE ID: "+te_id+" not found in IDs from GFF: "+locations+"\nplease make sure each ID in: "+taxonomy+" is found in:"+locations+"\n")
                    
                    if te_family not in consensus_seq_names:
                        sys.exit("TE Family: "+te_family+" not found in sequence names from: "+consensus+"\nplease make sure each family in: "+taxonomy+" is found in: "+consensus+"\n")
                

    #check coverage fasta
    if coverage_fasta is not None:
        mccutils.log("setup","checking coverage fasta: "+coverage_fasta)
        try:
            with open(coverage_fasta,"r") as fa:
                for record in SeqIO.parse(fa, "fasta"):
                    pass
        except Exception as e:
            print(e)
            sys.exit(coverage_fasta+" appears to be a malformed FastA file..exiting...\n")

    #check augment fasta
    if augment_fasta is not None:
        mccutils.log("setup","checking augment fasta: "+augment_fasta)
        try:
            with open(augment_fasta,"r") as fa:
                for record in SeqIO.parse(fa, "fasta"):
                    pass
        except Exception as e:
            print(e)
            sys.exit(augment_fasta+" appears to be a malformed FastA file..exiting...\n")       






def install(clean=False, debug=False):

    mcc_path = os.path.dirname(os.path.abspath(__file__))
    install_path = mcc_path+"/install/"
    install_config = install_path+"/config.json"
    log_dir = install_path+"/log/"
    conda_env_dir = install_path+"/envs/conda"
    data = {}
    data['paths'] = {
        'mcc_path': mcc_path,
        'install' : install_path,
        'log_dir': log_dir
    }
    
    data['URLs'] = config_install.URL
    data['MD5s'] = config_install.MD5
    data['ENVs'] = config_install.ENV
    data['output'] = config_install.OUTPUT

    for method in data['ENVs'].keys():
        data['ENVs'][method] = data['ENVs'][method].replace(config_install.ENV_PATH, install_path+"envs/")
    
    for method in data['output'].keys():
        data['output'][method] = data['output'][method].replace(config_install.INSTALL_PATH, install_path)

    with open(install_config,"w") as c:
        json.dump(data, c, indent=4)

    if os.path.exists(install_path+"install.log"):
        os.remove(install_path+"install.log")


    # removes installed tools and conda environments
    if clean:
        mccutils.log("install","Removing conda envs from: "+conda_env_dir)
        mccutils.log("install","Removing installed tools from: "+install_path+"tools")
        mccutils.remove(conda_env_dir)
        mccutils.remove(install_path+"/tools")

    mccutils.mkdir(conda_env_dir)
    os.chdir(install_path)
    mccutils.mkdir(log_dir)

    for env in config.ALL_METHODS:
        if env not in config.NO_INSTALL_METHODS:
            mccutils.log("install","Installing conda environment for: "+env)
            command = ["snakemake","--use-conda", "--conda-frontend","mamba", "--conda-prefix", conda_env_dir, "--configfile", install_config, "--cores", "1", "--nolock", "--conda-create-envs-only", data['output'][env]]

            if not debug:
                command.append("--quiet")
            mccutils.run_command(command)

            mccutils.log("install","Installing scripts for:"+env)
            command = ["snakemake","--use-conda", "--conda-prefix", conda_env_dir, "--configfile", install_config, "--cores", "1", "--nolock", data['output'][env]]
            if not debug:
                command.append("--quiet")
            mccutils.run_command(command)
    
    mccutils.log("install", "Installing conda environment for processing steps")
    command = ["snakemake","--use-conda", "--conda-frontend","mamba", "--conda-prefix", conda_env_dir, "--configfile", install_config, "--cores", "1", "--nolock", "--conda-create-envs-only", data['output']['processing']]

    if not debug:
        command.append("--quiet")
    mccutils.run_command(command)


def make_run_config(args, sample_name, ref_name, full_command, current_directory):
    run_id = random.randint(1000000,9999999)
    mccutils.mkdir(args.out+"/snakemake")
    mccutils.mkdir(args.out+"/snakemake/config")
    run_config = args.out+"/snakemake/config/config_"+str(run_id)+".json"
    input_dir = args.out
    reference_dir = args.out+"/"+ref_name+"/"
    sample_dir = args.out+"/"+sample_name+"/"
    results_dir = args.out+"/"+sample_name+"/results/"

    mcc_path = os.path.dirname(os.path.abspath(__file__))

    # get git commit hash to provide in summary report
    git_commit = "?"
    try:
        os.chdir(mcc_path)
        git_commit_file = args.out+"/git-commit.txt"
        mccutils.run_command_stdout(["git","rev-parse","HEAD"], git_commit_file)
        with open(git_commit_file,"r") as inf:
            for line in inf:
                git_commit = line.replace("\n","")
        
        mccutils.remove(git_commit_file)

    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("Could not locate git commit hash...using '?' ", file=sys.stderr)
        git_commit = "?"

    mccutils.log("SETUP","McClintock Version: "+git_commit)

    out_files_to_make = []
    out_files = config.OUT_PATHS
    for key in out_files.keys():
        out_files[key] = out_files[key].replace(config.INPUT_DIR, input_dir)
        out_files[key] = out_files[key].replace(config.REF_DIR, reference_dir)
        out_files[key] = out_files[key].replace(config.SAM_DIR, sample_dir)
        out_files[key] = out_files[key].replace(config.RESULTS_DIR, results_dir)
        out_files[key] = out_files[key].replace(config.SAMPLE_NAME, sample_name)   
    
    for method in args.methods:
        out_files_to_make.append(out_files[method])

    now = datetime.now()
    now_str = now.strftime("%Y%m%d.%H%M%S")
    log_dir = args.out+"/logs/"+now_str+"."+str(run_id)+"/"
    mccutils.mkdir(log_dir)


    chromosomes = []
    for record in SeqIO.parse(args.reference, "fasta"):
        chrom = str(record.id)
        chrom = mccutils.replace_special_chars(chrom)
        chromosomes.append(chrom)


    data = {}
    data['args'] = {
        'proc': str(args.proc),
        'out': sample_dir,
        'log_dir': log_dir,
        'augment_fasta': str(args.augment),
        'mcc_path': mcc_path,
        'commit': git_commit,
        'sample_name': sample_name,
        'ref_name': ref_name,
        'run_id' : str(run_id),
        'methods' : ",".join(args.methods),
        'out_files': ",".join(out_files_to_make),
        'save_comments' : str(args.comments),
        'max_threads_per_rule' : max(1, calculate_max_threads(args.proc, args.methods, config.MULTI_THREAD_METHODS, slow=args.slow)),
        'full_command' : full_command,
        'call_directory': current_directory,
        'time': now.strftime("%Y-%m-%d %H:%M:%S"),
        "chromosomes" : ",".join(chromosomes)
    }

    # input paths for files
    data["in"] = {
        'reference' : str(args.reference),
        'consensus' : str(args.consensus),
        'fq1': str(args.first),
        'fq2': str(args.second),
        'locations': str(args.locations),
        'taxonomy': str(args.taxonomy),
        'coverage_fasta': str(args.coverage_fasta),
    }

    # where mcc copies will be stored
    
    data["mcc"] = config.INTERMEDIATE_PATHS
    for key in data["mcc"].keys():
        data["mcc"][key] = data["mcc"][key].replace(config.INPUT_DIR, input_dir)
        data["mcc"][key] = data["mcc"][key].replace(config.REF_DIR, reference_dir)
        data["mcc"][key] = data["mcc"][key].replace(config.SAM_DIR, sample_dir)
        data["mcc"][key] = data["mcc"][key].replace(config.REF_NAME, ref_name)
        data["mcc"][key] = data["mcc"][key].replace(config.SAMPLE_NAME, sample_name)
    
    data["out"] = config.OUT_PATHS
    for key in data["mcc"].keys():
        data["mcc"][key] = data["mcc"][key].replace(config.INPUT_DIR, input_dir)
        data["mcc"][key] = data["mcc"][key].replace(config.REF_DIR, reference_dir)
        data["mcc"][key] = data["mcc"][key].replace(config.SAM_DIR, sample_dir)
        data["mcc"][key] = data["mcc"][key].replace(config.REF_NAME, ref_name)
        data["mcc"][key] = data["mcc"][key].replace(config.SAMPLE_NAME, sample_name)
    

    env_path = os.path.dirname(os.path.abspath(__file__))+"/install/envs/"
    data["envs"] = config_install.ENV
    for key in data["envs"].keys():
        data['envs'][key] = data['envs'][key].replace(config_install.ENV_PATH, env_path)


    with open(run_config,"w") as conf:
        json.dump(data, conf, indent=4)
    
    return run_id


def run_workflow(args, sample_name, ref_name, run_id, debug=False, annotations_only=False):
    log = args.out+"/mcclintock."+str(run_id)+".log"

    input_dir = args.out
    reference_dir = args.out+"/"+ref_name+"/"
    sample_dir = args.out+"/"+sample_name+"/"
    results_dir = args.out+"/"+sample_name+"/results/"

    out_files = config.OUT_PATHS
    for key in out_files.keys():
        out_files[key] = out_files[key].replace(config.INPUT_DIR, input_dir)
        out_files[key] = out_files[key].replace(config.REF_DIR, reference_dir)
        out_files[key] = out_files[key].replace(config.SAM_DIR, sample_dir)
        out_files[key] = out_files[key].replace(config.RESULTS_DIR, results_dir)
        out_files[key] = out_files[key].replace(config.SAMPLE_NAME, sample_name)    

    path=os.path.dirname(os.path.abspath(__file__))
    mccutils.mkdir(args.out+"/snakemake")
    snakemake_path = args.out+"/snakemake/"+str(run_id)
    mccutils.mkdir(snakemake_path)
    mccutils.run_command(["cp", path+"/Snakefile", snakemake_path])
    os.chdir(snakemake_path)
    command = ["snakemake","--use-conda", "--conda-prefix", path+"/install/envs/conda"]
    if not debug:
        command.append("--quiet")
    else:
        command.append("--reason")
    
    command += ["--configfile", args.out+"/snakemake/config/config_"+str(run_id)+".json"]
    command += ["--cores", str(args.proc)]

    mccutils.mkdir(sample_dir)
    mccutils.mkdir(sample_dir+"tmp")

    if args.clean:
        clean_command = command + ["--delete-all-output"]
        mccutils.run_command(clean_command)
        mccutils.remove(args.out+"/input")


    if not annotations_only:
        for method in args.methods:
            command.append(out_files[method])
        
        command.append(sample_dir+"results/summary/data/run/summary_report.txt")
    else:
        command.append(reference_dir+"reference_te_locations/inrefTEs.gff")
        command.append(reference_dir+"te_taxonomy/taxonomy.tsv")

    # print(" ".join(command))
    try:
        mccutils.run_command(command)
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("McClintock Pipeline Failed... please open an issue at https://github.com/bergmanlab/mcclintock/issues if you are having trouble using McClintock", file=sys.stderr)
        sys.exit(1)
    mccutils.remove(sample_dir+"tmp")


def calculate_max_threads(avail_procs, methods_used, multithread_methods, slow=False):
    max_threads = avail_procs

    if not slow:
        multi_methods_used = 0
        for method in methods_used:
            if method in multithread_methods:
                multi_methods_used += 1

        if multi_methods_used > 1:
            is_even = False
            if max_threads%2 == 0:
                is_even = True
            
            max_threads = max_threads//2
            if is_even:
                max_threads = max_threads-1

    return max_threads        


if __name__ == "__main__":                
    main()