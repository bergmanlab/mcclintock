import os
import sys
import subprocess
import traceback
import importlib.util as il
spec = il.spec_from_file_location("config", snakemake.params.config)
config = il.module_from_spec(spec)
sys.modules[spec.name] = config
spec.loader.exec_module(config)
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    reference_te_gff = snakemake.input.reference_tes
    bam = snakemake.input.bam

    out_dir = snakemake.params.out_dir
    script_dir = snakemake.params.script_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    status_log = snakemake.params.status_log
    threads = snakemake.threads

    out = snakemake.output.out

    mccutils.log("jitterbug","Running jitterbug", log=log)

    try:
        out_gff, config_file = run_jitterbug(
                                script_dir, 
                                bam, 
                                reference_te_gff, 
                                sample_name, out_dir, 
                                minmapq=config.RUN['MINMAPQ'], 
                                min_cluster_size=config.RUN['MIN_CLUSTER_SIZE'], 
                                threads=threads, 
                                log=log
        )

        config_file = make_config(
                            config_file, 
                            out_dir,
                            cluster_size=config.FILTER["CLUSTER_SIZE"],
                            span=config.FILTER['SPAN'],
                            int_size=config.FILTER['INT_SIZE'],
                            softclipped=config.FILTER['SOFTCLIPPED'],
                            pick_consistent=config.FILTER['PICK_CONSISTENT']
        )
        filter_jitterbug(script_dir, out_gff, config_file, sample_name, out, log=log)
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        with open(log,"a") as l:
            print(track, file=l)
        mccutils.log("Jitterbug","Jitterbug run failed")
        with open(status_log,"w") as l:
            l.write("FAILED\n")

        mccutils.run_command(["touch", out])


def run_jitterbug(script_dir, bam, ref_te_gff, sample_name, out_dir, minmapq=15, min_cluster_size=2, threads=1, log=None):
    command = [
        script_dir+"jitterbug.py", 
            "--mem", 
            "--numCPUs", str(threads), 
            "--output_prefix", out_dir+"/"+sample_name, 
            "-q", str(minmapq),
            "-c", str(min_cluster_size),
            bam,
            ref_te_gff
    ]
    mccutils.run_command(command, log=log)

    return out_dir+"/"+sample_name+".TE_insertions_paired_clusters.gff3", out_dir+"/"+sample_name+".filter_config.txt"


def make_config(config_file, out_dir, cluster_size=None, span=None, int_size=None, softclipped=None, pick_consistent=None):
    if cluster_size == None:
        out_cluster = [None, None]
    else:
        out_cluster = cluster_size

    if span == None:
        out_span = [None, None]
    else:
        out_span = span
    
    if int_size == None:
        out_int_size = [None, None]
    else:
        out_int_size = int_size
    
    if softclipped == None:
        out_softclipped = [None, None]
    else:
        out_softclipped = softclipped
    
    if pick_consistent == None:
        out_pick_consistent = [None, None]
    else:
        out_pick_consistent = pick_consistent
    
    with open(config_file, "r") as config:
        for line in config:
            line = line.replace("\n","")
            split_line = line.split("\t")
            if "cluster_size" in split_line[0]:
                if out_cluster[0] == None:
                    out_cluster[0] = split_line[1]
                
                if out_cluster[1] == None:
                    out_cluster[1] = split_line[2]
            
            if "span" in split_line[0]:
                if out_span[0] == None:
                    out_span[0] = split_line[1]
                
                if out_span[1] == None:
                    out_span[1] = split_line[2]
            
            if "int_size" in split_line[0]:
                if out_int_size[0] == None:
                    out_int_size[0] = split_line[1]
                
                if out_int_size[1] == None:
                    out_int_size[1] = split_line[2]
            
            if "softclipped" in split_line[0]:
                if out_softclipped[0] == None:
                    out_softclipped[0] = split_line[1]
                
                if out_softclipped[1] == None:
                    out_softclipped[1] = split_line[2]
            
            if "pick_consistent" in split_line[0]:
                if out_pick_consistent[0] == None:
                    out_pick_consistent[0] = split_line[1]
                
                if out_pick_consistent[1] == None:
                    out_pick_consistent[1] = split_line[2]

    out_config = out_dir+"mcc.config.txt"
    with open(out_config, "w") as out:
        out.write("\t".join(["cluster_size", str(out_cluster[0]), str(out_cluster[1])])+"\n")
        out.write("\t".join(["span", str(out_span[0]), str(out_span[1])])+"\n")
        out.write("\t".join(["int_size", str(out_int_size[0]), str(out_int_size[1])])+"\n")
        out.write("\t".join(["softclipped", str(out_softclipped[0]), str(out_softclipped[1])])+"\n")
        out.write("\t".join(["pick_consistent", str(out_pick_consistent[0]), str(out_pick_consistent[1])])+"\n")
    
    return out_config

def filter_jitterbug(script_dir, jitterbug_gff, filter_config, sample_name, filtered_gff, log=None):
    command = [
        script_dir+"tools/jitterbug_filter_results_func.py", 
            "-g", jitterbug_gff, 
            "-c", filter_config,
            "-o", filtered_gff
    ]

    mccutils.run_command(command, log=log)

if __name__ == "__main__":                
    main()