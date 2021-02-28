#!/usr/bin/env python3

import argparse
import os
import sys
import json
import random
import gzip
import hashlib
from datetime import datetime
import traceback

try:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    import scripts.mccutils as mccutils
    import config._internal.config as internal_config
    from Bio import SeqIO
except ImportError as e:
    print(e)
    sys.exit("ERROR...unable to load required python modules\ntry reinstalling and/or activating mcclintock environment:\n\tconda env create -f install/envs/mcclintock.yml --name mcclintock\n\tconda activate mcclintock\n")





def main():
    full_command = " ".join(["python3"] + sys.argv)
    current_directory = os.getcwd()

    expected_configs = internal_config.CONFIGS
    args = parse_args(expected_configs)
    sys.path = [args.config] + sys.path
    import _internal.config as config
    import install.install as config_install

    mccutils.mkdir(args.out+"/logs")
    mccutils.mkdir(args.out+"/tmp")
    check_installed_modules(args.methods, config.NO_INSTALL_METHODS, config_install.MD5, os.path.dirname(os.path.abspath(__file__))+"/install/")
    check_input_files(args.reference, args.consensus, args.first, fq2=args.second, locations=args.locations, taxonomy=args.taxonomy, coverage_fasta=args.coverage_fasta, augment_fasta=args.augment, annotations_only=args.make_annotations)
    ref_name = mccutils.get_base_name(args.reference)
    run_id, out_files = make_run_config(args, args.sample_name, ref_name, full_command, current_directory, config, config_install, debug=args.debug)
    run_workflow(args, args.sample_name, ref_name, run_id, config, out_files, debug=args.debug, annotations_only=args.make_annotations)
    mccutils.remove(args.out+"/tmp")

def parse_args(expected_configs):
    parser = argparse.ArgumentParser(prog='McClintock', description="Meta-pipeline to identify transposable element insertions using next generation sequencing data")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required=('--install' not in sys.argv))
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required='--install' not in sys.argv)
    parser.add_argument("-1", "--first", type=str, help="The path of the first fastq file from paired end read sequencing or the fastq file from single read sequencing", required=(('--install' not in sys.argv) and ('--make_annotations' not in sys.argv)))
    

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
    parser.add_argument("--sample_name", type=str, help="The sample name to use for output files [default: fastq1 name]", required=False)
    parser.add_argument("--resume", action="store_true", help="This option will attempt to use existing intermediate files from a previous McClintock run", required=False)
    parser.add_argument("--install", action="store_true", help="This option will install the dependencies of mcclintock", required=False)
    parser.add_argument("--debug", action="store_true", help="This option will allow snakemake to print progress to stdout", required=False)
    parser.add_argument("--slow", action="store_true", help="This option runs without attempting to optimize thread usage to run rules concurrently. Each multithread rule will use the max processors designated by -p/--proc", required=False)
    parser.add_argument("--make_annotations", action="store_true", help="This option will only run the pipeline up to the creation of the repeat annotations", required=False)
    parser.add_argument("-k","--keep_intermediate", type=str, help="This option determines which intermediate files are preserved after McClintock completes [default: general][options: minimal, general, methods, <list,of,methods>, all]", required=False)
    parser.add_argument("--config", type=str, help="This option determines which config files to use for your mcclintock run [default: config in McClintock Repository]", required=False)

    args = parser.parse_args()

    if args.config is None:
        args.config = os.path.dirname(os.path.abspath(__file__)) + "/config/"
    else:
        args.config = os.path.abspath(args.config)+"/"

    for key in expected_configs.keys():
        for config_file in expected_configs[key]:
            if not os.path.exists(args.config+"/"+config_file):
                sys.exit("Error: can't find config file: "+args.config+"/"+config_file+"\n Check that --config is set correctly...exiting...\n")

    sys.path = [args.config] + sys.path
    import _internal.config as config
    import install.install as config_install

    if args.debug is None:
        args.debug = False

    #check -m
    # If only one fastq has been supplied assume this is single ended data and launch only ngs_te_mapper and RelocaTE
    if args.second is None and not args.install:
        valid_methods = config.SINGLE_END_METHODS #from config.py
    else:
        valid_methods = config.ALL_METHODS #from config.py
    
    # used to preserve trimgalore and mapped reads output if they are explicitly called by the user
    trimgalore_called = False
    map_reads_called = False

    if args.methods is None:
        args.methods = valid_methods
    
    else:
        args.methods = args.methods.split(",")
        if "trimgalore" in args.methods:
            trimgalore_called = True
        
        if "map_reads" in args.methods:
            map_reads_called = True
            
        for x,method in enumerate(args.methods):
            args.methods[x] = method.lower()
            if args.methods[x] not in valid_methods:
                sys.stderr.write(" ".join(["Method:",method, "not a valid method...", "Valid methods:"," ".join(valid_methods),"\n"]))
                sys.exit(1)

    if args.install:
        mccutils.log("install","installing dependencies")
        mccutils.log("install","WARNING: this could take awhile")
        install(args.methods, config_install, resume=args.resume, debug=args.debug)
        sys.exit(0)

    #check -r
    args.reference = mccutils.get_abs_path(args.reference)
    #check -c
    args.consensus = mccutils.get_abs_path(args.consensus)

    if args.make_annotations != True:
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
    
    # check sample name
    if args.sample_name is not None:
        if "/" in args.sample_name or args.sample_name == "tmp":
            sys.exit(args.sample_name+" is not a valid sample name...\n")
    else:
        if not args.make_annotations:
            args.sample_name = mccutils.get_base_name(args.first)
        else:
            args.sample_name = "tmp"

    keep_intermediate_options = ["minimal","general", "methods", "all"] + args.methods
    if args.keep_intermediate is None:
        args.keep_intermediate = ["general"]
    else:
        args.keep_intermediate = args.keep_intermediate.split(",")
        for option in args.keep_intermediate:
            if option not in keep_intermediate_options:
                sys.stderr.write("keep_intermediate option: "+option+" is not valid. Valid options: "+" ".join(keep_intermediate_options)+"\nExample:(--keep_intermediate general,methods)\n")
                sys.exit(1)

    if trimgalore_called:
        args.keep_intermediate.append("trimgalore")
    
    if map_reads_called:
        args.keep_intermediate.append("map_reads")

    return args

def check_input_files(ref, consensus, fq1, fq2=None, locations=None, taxonomy=None, coverage_fasta=None, augment_fasta=None, annotations_only=False):
    # check reference fasta
    format_fasta(ref)

    if not annotations_only:
        #check fq1
        check_fastq(fq1)
        
        #check fq2
        if fq2 is not None:
            check_fastq(fq2)
    
    # checking consensus
    consensus_seq_names = format_fasta(consensus)
    
    #check locations gff
    gff_ids = []
    if locations is not None:
        gff_ids = format_gff(locations)
    
    # check taxonomy
    if taxonomy is not None:
        format_taxonomy(taxonomy, gff_ids, consensus_seq_names, consensus, locations)

    #check coverage fasta
    if coverage_fasta is not None:
        format_fasta(coverage_fasta)

    #check augment fasta
    if augment_fasta is not None:
        format_fasta(augment_fasta)


def format_fasta(in_fasta):
    mccutils.log("setup","checking fasta: "+in_fasta)
    seq_names = []
    try:
        with open(in_fasta,"r") as infa:
            for record in SeqIO.parse(infa, "fasta"):
                seq_name = str(record.id)
                if "#" in seq_name:
                    org_seq_name = seq_name
                    seq_name = seq_name[:(seq_name.find("#"))]
                    mccutils.log("setup", in_fasta+": replacing "+org_seq_name+" with "+seq_name+" for compatibility with RepeatMasker")

                masked_seq_name = mccutils.replace_special_chars(seq_name)
                if seq_name != masked_seq_name:
                    mccutils.log("setup", in_fasta+": ERROR problematic symbol in feature name: "+seq_name+" ... reformat this feature name for compatibility with McClintock")
                    print("Problematic symbols:"," ".join(mccutils.INVALID_SYMBOLS))
                    sys.exit(1)
                
                if masked_seq_name not in seq_names:
                    seq_names.append(masked_seq_name)
                else:
                    sys.exit(in_fasta+": Duplicate sequence name:"+masked_seq_name+"...exiting...\n")

    except Exception as e:
        print(e)
        sys.exit(in_fasta+" appears to be a malformed FastA file..exiting...\n")
    
    if len(seq_names) < 1:
        sys.exit(in_fasta+" contains no sequences... exiting...\n")

    return seq_names

def check_fastq(fastq):
    mccutils.log("setup","checking fastq: "+fastq)
    #check fq1
    if ".fastq" not in fastq and ".fq" not in fastq:
        sys.exit(fastq+" is not a (.fastq/.fq) file, exiting...\n")
    
    if not os.path.isfile(fastq):
        sys.exit(fastq+" does not exist... exiting...\n")
    
    if (mccutils.is_empty_file(fastq)):
        sys.exit(fastq+" is empty... exiting...\n")

def format_gff(ingff):
    mccutils.log("setup","checking locations gff: "+ingff)
    gff_ids = []
    with open(ingff,"r") as gff:
        for line in gff:
            if "#" not in line[0]:
                split_line = line.split("\t")
                if len(split_line) < 9:
                    sys.exit(ingff+" appears to be a malformed GFF file..exiting...\n")
                else:
                    feats = split_line[8]
                    split_feats = feats.split(";")
                    gff_id = ""
                    for feat in split_feats:
                        if feat[:3] == "ID=":
                            gff_id = feat.split("=")[1].replace("\n","")
                            masked_gff_id = mccutils.replace_special_chars(gff_id)
                            if gff_id != masked_gff_id:
                                mccutils.log("setup", ingff+": ERROR problematic symbol in feature name: "+gff_id+" ... reformat this feature name for compatibility with McClintock")
                                print("Problematic symbols:"," ".join(mccutils.INVALID_SYMBOLS))
                                sys.exit(1)

                            if masked_gff_id not in gff_ids:
                                gff_ids.append(masked_gff_id)
                            else:
                                sys.exit("ID: "+masked_gff_id+" is not unique. please ensure each feature has a unique ID\n")
                    if masked_gff_id == "":
                        sys.exit("GFF line: "+line+" is missing an ID attribute (ex. ID=chr1_TY1s1)\n")
    
    return gff_ids

def format_taxonomy(in_taxonomy, gff_ids, consensus_ids, consensus_fasta, locations_gff):
    mccutils.log("setup","checking taxonomy TSV: "+in_taxonomy)
    with open(in_taxonomy, "r") as tsv:
        for line in tsv:
            split_line = line.split("\t")
            if len(split_line) != 2:
                sys.exit(in_taxonomy+" does not have two columns. Should be tab-separated file with feature ID and TE family as columns\n")
            else:
                te_id = split_line[0]
                masked_te_id = mccutils.replace_special_chars(te_id)
                if masked_te_id != te_id:
                    mccutils.log("setup", in_taxonomy+": ERROR problematic symbol in feature name: "+te_id+" ... reformat this feature name for compatibility with McClintock")
                    print("Problematic symbols:"," ".join(mccutils.INVALID_SYMBOLS))
                    sys.exit(1)

                te_family = split_line[1].replace("\n","")
                if "#" in te_family:
                    org_te_family = te_family
                    te_family = te_family[:(te_family.find("#"))]
                    mccutils.log("setup", in_taxonomy+": replacing "+org_te_family+" with "+te_family+" for compatibility with RepeatMasker")

                masked_te_family = mccutils.replace_special_chars(te_family)
                if masked_te_family != te_family:
                    mccutils.log("setup", in_taxonomy+": ERROR problematic symbol in feature name: "+te_family+" ... reformat this feature name for compatibility with McClintock")
                    print("Problematic symbols:"," ".join(mccutils.INVALID_SYMBOLS))
                    sys.exit(1)

                if masked_te_id not in gff_ids:
                    sys.exit("TE ID: "+masked_te_id+" not found in IDs from GFF: "+locations_gff+"\nplease make sure each ID in: "+in_taxonomy+" is found in:"+locations_gff+"\n")
                
                if masked_te_family not in consensus_ids:
                    sys.exit("TE Family: "+masked_te_family+" not found in sequence names from: "+consensus_fasta+"\nplease make sure each family in: "+in_taxonomy+" is found in: "+consensus_fasta+"\n")

def get_conda_envs(conda_env_dir):
    existing_envs = {}
    if os.path.exists(conda_env_dir):
        for f in os.listdir(conda_env_dir):
            if ".yaml" in f:
                yaml = conda_env_dir+"/"+f
                name = ""
                with open(yaml,"r") as y:
                    for line in y:
                        if "name:" in line:
                            name = line.split(":")[1].replace("\n","")
                            name = name.replace(" ","")
                            existing_envs[name] = yaml
    
    return existing_envs

def install(methods, config_install, resume=False, debug=False):

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


    # finding existing conda yamls
    existing_envs = get_conda_envs(conda_env_dir)


    mccutils.mkdir(conda_env_dir)
    os.chdir(install_path)
    mccutils.mkdir(log_dir)

    # temp requires te-locate scripts to make taxonomy file
    if "temp" in methods and "te-locate" not in methods:
        methods.append("te-locate")

    for env in methods:
        if not resume:
            # remove existing envs
            if env in existing_envs.keys():
                mccutils.log("install","Removing existing conda env for: "+env)
                mccutils.remove(existing_envs[env])
                mccutils.remove(existing_envs[env].replace(".yaml",""))

            # remove existing src code
            if os.path.exists(install_path+"/tools/"+env):
                mccutils.log("install","Removing existing installation of: "+env)
                print(install_path+"/tools/"+env)
                mccutils.remove(install_path+"/tools/"+env)
        
        # reinstall src code
        mccutils.log("install","Installing scripts for:"+env)
        command = ["snakemake","--use-conda", "--conda-prefix", conda_env_dir, "--configfile", install_config, "--cores", "1", "--nolock", data['output'][env]]
        if not debug:
            command.append("--quiet")
        mccutils.run_command(command)

def get_installed_versions(tool_dir):
    versions = {}
    for d in os.listdir(tool_dir):
        version = ""
        version_file = tool_dir+"/"+d+"/version.log"
        if os.path.exists(version_file):
            with open(version_file,"r") as v_file:
                for line in v_file:
                    version = line
        
        versions[d] = version
    
    return versions


def check_installed_modules(methods, no_install_methods, method_md5s, install_dir):
    # finding existing conda yamls
    conda_env_dir = install_dir+"/envs/conda/"
    existing_envs = get_conda_envs(conda_env_dir)
    installed_version = get_installed_versions(install_dir+"/tools/")

    for method in methods:
        if method not in existing_envs.keys():
            print("ERROR: missing conda env for method:", method, file=sys.stderr)
            print("please install methods using: python mcclintock.py --install", file=sys.stderr)
            sys.exit(1)
        
        if method not in no_install_methods:
            if installed_version[method] != method_md5s[method]:
                print("ERROR: installed version of", method,"(", installed_version[method], ") does not match the expected: ", method_md5s[method])
                print("please install methods using: python mcclintock.py --install", file=sys.stderr)
                sys.exit(1)


# stores config file info and hashes so that future runs can determine if config files have changed
def setup_config_info(config_dir, method_to_config, config_rules):
    out_dict = {}
    out_dict["path"] = config_dir
    for method in method_to_config.keys():
        out_dict[method] = {}
        out_dict[method]["files"] = []
        out_dict[method]["MD5s"] = []
        out_dict[method]['rules'] = config_rules[method]
        for f in method_to_config[method]:
            config_file = config_dir+f
            with open(config_file,"rb") as cf:
                data = cf.read()
                md5 = hashlib.md5(data).hexdigest()

            out_dict[method]["files"].append(config_file)
            out_dict[method]["MD5s"].append(md5)
    
    return out_dict

def make_run_config(args, sample_name, ref_name, full_command, current_directory, config, config_install, debug=False):
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
        passed = mccutils.run_command_stdout(["git","rev-parse","HEAD"], git_commit_file, fatal=False)
        if not passed:
            raise Exception("Could not locate git commit hash")
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

    method_paths = config.OUT_DIRS
    for method in method_paths.keys():
        method_paths[method] = method_paths[method].replace(config.RESULTS_DIR, results_dir)
        method_paths[method] = method_paths[method].replace(config.SAM_DIR, sample_dir)

    out_files_to_make = []
    out_files = config.OUT_PATHS
    for key in out_files.keys():
        out_files[key] = out_files[key].replace(config.METHOD_DIR, method_paths[key])
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
        'run_config': run_config,
        'methods' : ",".join(args.methods),
        'out_files': ",".join(out_files_to_make),
        'save_comments' : str(args.comments),
        'max_threads_per_rule' : max(1, calculate_max_threads(args.proc, args.methods, config.MULTI_THREAD_METHODS, slow=args.slow)),
        'full_command' : full_command,
        'call_directory': current_directory,
        'time': now.strftime("%Y-%m-%d %H:%M:%S"),
        "chromosomes" : ",".join(chromosomes),
        "debug": str(debug)
    }

    data["config"] = setup_config_info(args.config, config.CONFIGS, config.CONFIG_RULES)

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
    
    data["out"] = out_files

    data['outdir'] = method_paths
    
    data["essential"] = config.ESSENTIAL_PATHS
    for key in data["essential"].keys():
        for x,val in enumerate(data["essential"][key]):
            data["essential"][key][x] = val.replace(config.METHOD_DIR, method_paths[key])
            data["essential"][key][x] = data["essential"][key][x].replace(config.SAMPLE_NAME, sample_name)
            data["essential"][key][x] = data["essential"][key][x].replace(config.REF_NAME, ref_name)

    env_path = os.path.dirname(os.path.abspath(__file__))+"/install/envs/"
    data["envs"] = config_install.ENV
    for key in data["envs"].keys():
        data['envs'][key] = data['envs'][key].replace(config_install.ENV_PATH, env_path)


    with open(run_config,"w") as conf:
        json.dump(data, conf, indent=4)
    
    return run_id, out_files


def run_workflow(args, sample_name, ref_name, run_id, config, out_files, debug=False, annotations_only=False):
    log = args.out+"/mcclintock."+str(run_id)+".log"

    input_dir = args.out
    reference_dir = args.out+"/"+ref_name+"/"
    sample_dir = args.out+"/"+sample_name+"/"
    results_dir = args.out+"/"+sample_name+"/results/" 

    path=os.path.dirname(os.path.abspath(__file__))
    mccutils.mkdir(args.out+"/snakemake")
    snakemake_path = args.out+"/snakemake/"+str(run_id)
    mccutils.mkdir(snakemake_path)
    mccutils.run_command(["cp", path+"/Snakefile", snakemake_path])
    os.chdir(snakemake_path)
    command = ["snakemake","--use-conda", "--conda-prefix", path+"/install/envs/conda"]

    config_json = args.out+"/snakemake/config/config_"+str(run_id)+".json"
    command += ["--configfile", config_json]
    command += ["--cores", str(args.proc)]

    if not args.resume:
        if os.path.exists(reference_dir) and len(os.listdir(reference_dir)) > 0:
            mccutils.remove(config_json)
            sys.exit("ERROR: output directory:"+reference_dir+" is not empty. If wanting to resume a previous run, use --resume, otherwise please delete this directory or change your -o/--output\n")
        if os.path.exists(sample_dir) and len(os.listdir(sample_dir)) > 0:
            mccutils.remove(config_json)
            sys.exit("ERROR: output directory:"+sample_dir+" is not empty. If wanting to resume a previous run, use --resume, otherwise please delete this directory or change your -o/--output or --sample_name\n")

    
    # check that previous runs are compatible
    else:
        mccutils.log("setup","Checking config files to ensure previous intermediate files are compatible with this run")
        config_found = False
        previous_config_md5s = {}
        for prev_config in os.listdir(input_dir+"/snakemake/config/"):
            if prev_config != "config_"+str(run_id)+".json":
                config_found = True
                config_compatible = config_compatibility(input_dir+"/snakemake/config/config_"+str(run_id)+".json", args.out+"/snakemake/config/"+prev_config)
                previous_config_md5s = get_recent_config_md5s(args.out+"/snakemake/config/"+prev_config, previous_config_md5s)
                if not config_compatible:
                    mccutils.remove(config_json)
                    sys.exit(1)
        
        if not config_found:
            mccutils.remove(config_json)
            sys.exit("ERROR: Unable to resume run. No config files from previous runs found in:"+input_dir+"/snakemake/config/ Remove --resume for clean run\n")

        rules_to_rerun = get_rules_to_rerun(input_dir+"/snakemake/config/config_"+str(run_id)+".json", previous_config_md5s)
        command.append("-R")
        if len(rules_to_rerun) > 0:
            for rule in rules_to_rerun:
                command.append(rule)

    if not debug:
        command.append("--quiet")
    else:
        command.append("--reason")
        command.append("--verbose")

    if not annotations_only:
        for method in args.methods:
            command.append(out_files[method])
        
        command.append(sample_dir+"results/summary/data/run/summary_report.txt")
    else:
        command.append(reference_dir+"reference_te_locations/inrefTEs.gff")
        command.append(reference_dir+"te_taxonomy/taxonomy.tsv")

    # print(" ".join(command))
    try:
        sys.stdout.flush()
        mccutils.mkdir(sample_dir)
        mccutils.mkdir(sample_dir+"tmp")
        if debug:
            print(" ".join(command))
        mccutils.run_command(command)
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("McClintock Pipeline Failed... please open an issue at https://github.com/bergmanlab/mcclintock/issues if you are having trouble using McClintock", file=sys.stderr)
        sys.exit(1)
    mccutils.remove(sample_dir+"tmp")
    remove_intermediate_files(args.keep_intermediate, config_json, args.methods, ref_name, sample_name, args.out)

def get_recent_config_md5s(prev_config, config_md5s):
    with open(prev_config) as f:
        prev_config_data = json.load(f)
    
    run_methods = prev_config_data['args']['methods']
    start_time = datetime.strptime(prev_config_data['args']['time'], '%Y-%m-%d %H:%M:%S')

    for config in prev_config_data['config'].keys():
        if config in run_methods:
            if config not in config_md5s.keys():
                config_md5s[config] = [prev_config_data['config'][config]['MD5s'], start_time]
            else:
                if config_md5s[config][1] < start_time:
                    config_md5s[config] = [prev_config_data['config'][config]['MD5s'], start_time]
    
    return config_md5s

def get_rules_to_rerun(run_config, prev_md5s):
    rules_to_rerun = []
    with open(run_config) as f:
        run_config_data = json.load(f)
    
    for config in run_config_data['config'].keys():
        if config in prev_md5s.keys():
            for x,md5 in enumerate(prev_md5s[config][0]):
                if md5 != run_config_data['config'][config]["MD5s"][x]:
                    rules_to_rerun.append(run_config_data['config'][config]['rules'][x])

    return rules_to_rerun

def config_compatibility(run_config, prev_config):
    with open(run_config) as f:
        run_config_data = json.load(f)
    
    with open(prev_config) as f:
        prev_config_data = json.load(f)
    
    # check McClintock commit compatibility
    if run_config_data['args']['commit'] != prev_config_data['args']['commit']:
        sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
        sys.stderr.write("                  Previous mcclintock run intermediate files are incompatible due to differences in McClintock commit version\n")
        sys.stderr.write("                  Current McClintock Version: "+run_config_data['args']['commit']+"\n")
        sys.stderr.write("                  Previous McClintock run Version: "+prev_config_data['args']['commit']+"\n")
        return False

    # check ref compatibility
    if os.path.exists(run_config_data["mcc"]["reference"]):
        if run_config_data["mcc"]["reference"] == prev_config_data["mcc"]["reference"]:
            if run_config_data["args"]["augment_fasta"] != prev_config_data["args"]["augment_fasta"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous mcclintock run intermediate files are incompatible due to differences in --augment\n")
                return False
            
            if run_config_data["in"]["reference"] != prev_config_data["in"]["reference"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in -r/--reference\n")
                return False
    
    # check consensus compatibility
    if os.path.exists(run_config_data["mcc"]["consensus"]):
        if run_config_data["mcc"]["consensus"] == prev_config_data["mcc"]["consensus"] and run_config_data["in"]["consensus"] != prev_config_data["in"]["consensus"]:
            sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
            sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in -c/--consensus\n")
            return False
    
    # check fastq compatibility
    if "--make_annotations" not in prev_config_data['args']['full_command']:
        if os.path.exists(run_config_data["mcc"]["fq1"]):
            if run_config_data["in"]["fq1"] != prev_config_data["in"]["fq1"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in -1/--first\n")
                return False

            if run_config_data["in"]["fq2"] != prev_config_data["in"]["fq2"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in -2/--second\n")
                return False

            if ("trimgalore" in run_config_data["args"]["methods"] and "trimgalore" not in prev_config_data["args"]["methods"]) or ("trimgalore" in prev_config_data["args"]["methods"] and "trimgalore" not in run_config_data["args"]["methods"]):
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in running Trimgalore\n")
                return False
    
    # check locations gff compatibility
    if os.path.exists(run_config_data["mcc"]["locations"]):
        if run_config_data["mcc"]["locations"] == prev_config_data["mcc"]["locations"]:
            if run_config_data["args"]["augment_fasta"] != prev_config_data["args"]["augment_fasta"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous mcclintock run intermediate files are incompatible due to differences in --augment\n")
                return False

            if run_config_data["in"]["locations"] != prev_config_data["in"]["locations"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in -g/--locations\n")
                return False
            
    # check taxonomy compatibility
    if os.path.exists(run_config_data["mcc"]["taxonomy"]):
        if run_config_data["mcc"]["taxonomy"] == prev_config_data["mcc"]["taxonomy"]:
            if run_config_data["args"]["augment_fasta"] != prev_config_data["args"]["augment_fasta"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous mcclintock run intermediate files are incompatible due to differences in --augment\n")
                return False

            if run_config_data["in"]["taxonomy"] != prev_config_data["in"]["taxonomy"]:
                sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
                sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in -t/--taxonomy\n")
                return False

    # check coverage fasta compatibility
    if os.path.exists(run_config_data["mcc"]["coverage_fasta"]):
        if run_config_data["mcc"]["coverage_fasta"] == prev_config_data["mcc"]["coverage_fasta"] and run_config_data["in"]["coverage_fasta"] != prev_config_data["in"]["coverage_fasta"]:
            sys.stderr.write("(--resume) ERROR: Unable to resume McClintock run\n")
            sys.stderr.write("                  Previous McClintock run intermediate files are incompatible due to differences in -s/--coverage_fasta\n")
            return False


    return True


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


def remove_intermediate_files(options, run_config_file, methods, ref_name, sample_name, outdir):
    if "all" in options:
        return

    with open(run_config_file) as f:
        run_config_data = json.load(f)

    keep_paths = []
    if "methods" not in options:
        for method in methods:
            method_out = "/".join(run_config_data['out'][method].split("/")[:-1])+"/"
            if method not in options:
                essential_paths = run_config_data['essential'][method]
                if os.path.exists(method_out):
                    # delete all files not marked as essential
                    for root, subdirs, files in os.walk(method_out, topdown=False):
                        for f in files:
                            file_path = os.path.join(root, f)
                            is_essential = False
                            for essential_path in essential_paths:
                                if (os.path.isdir(essential_path) and essential_path in file_path) or (essential_path == file_path):
                                    is_essential = True
                            
                            if not is_essential:
                                mccutils.remove(file_path)
                    
                    # remove empty directories
                    for root, subdirs, files in os.walk(method_out, topdown=False):
                        for d in subdirs:
                            dir_path = os.path.join(root, d)
                            if len(os.listdir(dir_path)) < 1:
                                mccutils.remove(dir_path)
            
            else:
                keep_paths.append(method_out)
    
    if "general" not in options:
        intermediate_dir = outdir+"/"+sample_name+"/intermediate/"
        for root, subdirs, files in os.walk(intermediate_dir, topdown=False):
            for f in files:
                file_path = os.path.join(root, f)
                keep = False
                for keep_path in keep_paths:
                    if keep_path in file_path:
                        keep = True
                
                if not keep:
                    mccutils.remove(file_path)
        
        # remove empty directories
        for root, subdirs, files in os.walk(intermediate_dir, topdown=False):
            for d in subdirs:
                dir_path = os.path.join(root, d)
                if len(os.listdir(dir_path)) < 1:
                    mccutils.remove(dir_path)

        if len(os.listdir(intermediate_dir)) < 1:
            mccutils.remove(intermediate_dir)







if __name__ == "__main__":                
    main()