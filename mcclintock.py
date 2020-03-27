#!/usr/bin/env python3

import argparse
import os
import modules.mccutils as mccutils
import config.config as config
import json
import random
import sys


def main():
    args = parse_args()
    mccutils.mkdir(args.out+"/logs")
    mccutils.mkdir(args.out+"/tmp")
    sample_name = mccutils.get_base_name(args.first, fastq=True)
    ref_name = mccutils.get_base_name(args.reference)
    run_id = make_run_config(args, sample_name, ref_name)
    run_workflow(args, sample_name, run_id)

def parse_args():
    parser = argparse.ArgumentParser(prog='McClintock', description="Meta-pipeline to identify transposable element insertions using next generation sequencing data")

    ## required ##
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required='--install' not in sys.argv)
    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs for the species in fasta format", required='--install' not in sys.argv)
    parser.add_argument("-1", "--first", type=str, help="The path of the first fastq file from paired end read sequencing or the fastq file from single read sequencing", required='--install' not in sys.argv)
    

    ## optional ##
    parser.add_argument("-2", "--second", type=str, help="The path of the second fastq file from a paired end read sequencing", required=False)
    parser.add_argument("-p", "--proc", type=int, help="The number of processors to use for parallel stages of the pipeline [default = 1]", required=False)
    parser.add_argument("-M", "--mem", type=int, help="The amount of memory available for the pipeline in GB. [default = 4]", required=False)
    parser.add_argument("-o", "--out", type=str, help="An output folder for the run. [default = '.']", required=False)
    parser.add_argument("-m", "--methods", type=str, help="A comma-delimited list containing the software you want the pipeline to use for analysis. e.g. '-m relocate,TEMP,ngs_te_mapper' will launch only those three methods", required=False)
    parser.add_argument("-g", "--locations", type=str, help="The locations of known TEs in the reference genome in GFF 3 format. This must include a unique ID attribute for every entry", required=False)
    parser.add_argument("-t", "--taxonomy", type=str, help="A tab delimited file with one entry per ID in the GFF file and two columns: the first containing the ID and the second containing the TE family it belongs to. The family should correspond to the names of the sequences in the consensus fasta file", required=False)
    parser.add_argument("-s", "--coverage_fasta", type=str, help="A fasta file that will be used for TE-based coverage analysis, if not supplied then the consensus sequences of the TEs will be used for the analysis", required=False)
    # parser.add_argument("-d", "--coverage", action="store_true", help="If this option is specified then McClintock will perform depth of coverage analysis for every TE. Note: Doing TE-based coverage analysis will result in longer running time. A fasta file can be provided here for coverage analysis. If no file is provided here, the consensus sequences of the TEs will be used for the analysis", required=False)
    # parser.add_argument("-D", "--coverage_only", action="store_true", help="If this option is specified then only depth of coverage analysis for TEs will be performed", required=False)
    parser.add_argument("-T", "--comments", action="store_true", help="If this option is specified then fastq comments (e.g. barcode) will be incorporated to SAM output. Warning: do not use this option if the input fastq files do not have comments", required=False)
    # parser.add_argument("-b", "--keep_bam", action="store_true", help="Retain the sorted and indexed BAM file of the paired end data aligned to the reference genome", required=False)
    # parser.add_argument("-i", "--remove_intermediate", action="store_true", help="If this option is specified then all sample specific intermediate files will be removed, leaving only the overall results. The default is to leave sample specific intermediate files", required=False)
    parser.add_argument("-C", "--include_consensus", action="store_true", help="This option will include the consensus TE sequences as extra chromosomes in the reference file (useful if the organism is known to have TEs that are not present in the reference strain)", required=False)
    parser.add_argument("-R", "--include_reference", action="store_true", help="This option will include the reference TE sequences as extra chromosomes in the reference file", required=False)
    parser.add_argument("--clean", action="store_true", help="This option will make sure mcclintock runs from scratch and doesn't reuse files already created", required=False)
    parser.add_argument("--install", action="store_true", help="This option will install the dependencies of mcclintock", required=False)

    args = parser.parse_args()

    if args.install:
        print("installing dependencies...")
        install(clean=args.clean)
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
    
    #check -M
    if args.mem is None:
        args.mem = 4

    #check -o
    if args.out is None:
        args.out = os.path.abspath(".")
    else:
        args.out = os.path.abspath(args.out)
        mccutils.mkdir(args.out)
    

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
    
    return args


def install(clean=False):

    mcc_path = os.path.dirname(os.path.abspath(__file__))
    install_path = mcc_path+"/install/"
    install_config = install_path+"/config.json"
    data = {}
    data['paths'] = {
        'mcc_path': mcc_path,
        'install' : install_path,
        'envs' : mcc_path+"/envs/"
    }
    with open(install_config,"w") as config:
        json.dump(data, config, indent=4)

    if os.path.exists(install_path+"install.log"):
        os.remove(install_path+"install.log")

    conda_env_dir = mcc_path+"/envs/conda"

    # removes installed tools and conda environments
    if clean:
        print("Removing conda envs from: "+conda_env_dir)
        print("Removing installed tools from: "+install_path+"tools")
        mccutils.run_command(["rm","-r",conda_env_dir, install_path+"/tools"])

    mccutils.mkdir(conda_env_dir)
    os.chdir(install_path)
    command = ["snakemake","--use-conda", "--conda-prefix", conda_env_dir, "--configfile", install_config, "-R", "install_all", "--cores", "1", "--nolock"]
    mccutils.run_command(command)

def make_run_config(args, sample_name, ref_name):
    run_id = random.randint(1000000,9999999)
    mccutils.mkdir(args.out+"/config")
    run_config = args.out+"/config/config_"+str(run_id)+".json"
    data = {}
    data['args'] = {
        'proc': str(args.proc),
        'mem': str(args.mem),
        'out': str(args.out),
        'include_consensus': str(args.include_consensus),
        'include_reference': str(args.include_reference),
        'mcc_path': os.path.dirname(os.path.abspath(__file__)),
        'sample_name': sample_name,
        'ref_name': ref_name,
        'run_id' : str(run_id),
        'methods' : ",".join(args.methods),
        'save_comments' : str(args.comments)
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
    input_dir = args.out+"/input/"
    data["mcc"] = {
        'mcc_files' : input_dir,
        'reference' : input_dir+ref_name+".fasta",
        'consensus' : input_dir+"consensusTEs.fasta",
        'fq1' : input_dir+sample_name+"_1.fastq",
        'fq2' : input_dir+sample_name+"_2.fastq",
        'locations' : input_dir+"inrefTEs.gff",
        'taxonomy' : input_dir+"taxonomy.tsv",
        'coverage_fasta' : input_dir+"coverageTEs.fasta",
        'masked_fasta' : input_dir+ref_name+".masked.fasta",
        'formatted_ref_tes' : input_dir+ref_name+".ref.TEs.gff",
        'formatted_taxonomy' : input_dir+ref_name+".TE.taxonomy.tsv",
        'formatted_consensus' : input_dir+"formattedConsensusTEs.fasta",
        'ref_te_fasta' : input_dir+ref_name+".ref.TEs.fasta",
        'augmented_reference' : input_dir+ref_name+".aug.fasta",
        'ref_tes_bed' : input_dir+ref_name+".ref.TEs.bed",
        'telocate_te_gff' : input_dir+ref_name+".ref.TEs_HL.gff",
        'telocate_sam' : input_dir+sample_name+".telocate.sam",
        'telocate_ref_fasta' : input_dir+ref_name+".aug.telocate.fasta",
        'sam' : input_dir+sample_name+".sam",
        'bam' : input_dir+sample_name+".bam",
        'ref_2bit' : input_dir+ref_name+".aug.fasta.2bit",
        'relocaTE_consensus' : input_dir+"formattedConsensusTEs.relocaTE.fasta",
        'relocaTE_ref_TEs' : input_dir+ref_name+".ref.TEs.relocaTE.gff",
        'popoolationTE_ref_fasta' : input_dir+ref_name+".masked.popoolationTE.fasta",
        'popoolationTE_taxonomy' : input_dir+"taxonomy.popoolationTE.tsv",
        'popoolationTE_gff' : input_dir+ref_name+".ref.TEs.popoolationTE.gff"
    }

    summary_dir = args.out+"/summary/"
    data['summary'] = {
        'flagstat' : input_dir+sample_name+".bam.flagstat",
        'median_insert_size' : input_dir+"median_insert.size"
    }


    with open(run_config,"w") as config:
        json.dump(data, config, indent=4)
    
    return run_id

def run_workflow(args, sample_name, run_id):
    log = args.out+"/mcclintock."+str(run_id)+".log"
    out_files = {
        'coverage': args.out+"/coverage/coverage.log",
        'ngs_te_mapper': args.out+"/ngs_te_mapper/ngs_te_mapper.log",
        'relocate': args.out+"/relocaTE/relocaTE.log",
        'temp': args.out+"/temp/temp.log",
        'retroseq': args.out+"/retroseq/retroseq.log",
        'popoolationte': args.out+"/popoolationTE/popoolationTE.log",
        'te-locate': args.out+"/te-locate/te-locate.log",
        'trimgalore': args.out+"/input/"+sample_name+"_1.fastq"
    }

    path=os.path.dirname(os.path.abspath(__file__))
    mccutils.mkdir(args.out+"/snakemake")
    snakemake_path = args.out+"/snakemake/"+str(run_id)
    mccutils.mkdir(snakemake_path)
    mccutils.run_command(["cp", path+"/Snakefile", snakemake_path])
    os.chdir(snakemake_path)
    command = ["snakemake","--use-conda", "--conda-prefix", path+"/envs/conda"]
    for method in args.methods:
        command.append(out_files[method])

    command += ["--configfile", args.out+"/config/config_"+str(run_id)+".json"]
    command += ["--cores", str(args.proc)]
    if args.clean:
        clean_command = command + ["--delete-all-output"]
        mccutils.run_command(clean_command)
        mccutils.run_command(["rm","-r",args.out+"/input"])

    print(" ".join(command))
    mccutils.run_command(command)


if __name__ == "__main__":                
    main()