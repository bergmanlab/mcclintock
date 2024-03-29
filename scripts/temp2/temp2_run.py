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
    fq1 = snakemake.input.fq1
    fq2 = snakemake.input.fq2
    bam = snakemake.input.bam
    reference = snakemake.input.reference
    twobit = snakemake.input.twobit
    consensus = snakemake.input.consensus
    ref_te_bed = snakemake.input.ref_te_bed
    taxonomy = snakemake.input.taxonomy
    median_insert_size_file = snakemake.input.median_insert_size
    log = snakemake.params.log

    with open(log,"a") as l:
        l.write("BAM: "+bam+"\n")
        l.write("2bit: "+twobit+"\n")
        l.write("consensus fasta: "+consensus+"\n")
        l.write("reference TE BED: "+ref_te_bed+"\n")
        l.write("Taxonomy TSV: "+taxonomy+"\n")

    threads = snakemake.threads
    out_dir = snakemake.params.out_dir
    script_dir = snakemake.params.script_dir
    sample_name = snakemake.params.sample_name
    status_log = snakemake.params.status_log

    # ensures intermediate files from previous runs are removed
    for f in os.listdir(out_dir):
        mccutils.remove(out_dir+"/"+f)
    
    mccutils.log("temp2","running TEMP2 Module")

    try:
        median_insert_size = get_median_insert_size(median_insert_size_file)
        run_temp2_insertion(fq1, fq2, bam, median_insert_size, reference, script_dir, consensus, ref_te_bed, threads, out_dir, config, log)
        run_temp2_absence(script_dir, bam, twobit, ref_te_bed, median_insert_size, threads, out_dir+"/absence", config, log)
        mccutils.run_command(["cp", out_dir+'/absence/'+sample_name+".absence.refined.bp.summary", out_dir], log=log)

        mccutils.check_file_exists(snakemake.output[0])
        mccutils.check_file_exists(snakemake.output[1])
        with open(status_log,"w") as l:
            l.write("COMPLETED\n")
        mccutils.log("temp2","TEMP2 run complete")
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        with open(log,"a") as l:
            print(track, file=l)
        mccutils.log("temp2","TEMP2 run failed")
        with open(status_log,"w") as l:
            l.write("FAILED\n")

        mccutils.run_command(["touch", snakemake.output[0]])
        mccutils.run_command(["touch", snakemake.output[1]])


def get_median_insert_size(infile):
    median_insert_size = 0
    with open(infile,"r") as inf:
        for line in inf:
            insert = line.split("=")[1]
            insert = insert.replace("\n","")
            median_insert_size = int(float(insert))
    
    return median_insert_size

def run_temp2_insertion(fq1, fq2, bam, insert_size, reference, scripts, consensus, te_bed, threads, out, config, log):
    mccutils.log("temp2","running TEMP2 non-reference insertion prediction", log=log)
    command = [
        "bash", scripts+"/TEMP2", "insertion", 
            "-l", fq1,
            "-r", fq2,
            "-i", bam,
            "-g", reference,
            "-R", consensus,
            "-t", te_bed,
            "-c", str(threads),
            "-f", str(insert_size),
            "-o", out
    ]
    
    for param in config.PARAMS['insertion'].keys():
        if config.PARAMS['insertion'][param] == True:
            command.append(param)
        elif config.PARAMS['insertion'][param] != False:
            command.append(param)
            command.append(str(config.PARAMS['insertion'][param]))

    mccutils.run_command(command, log=log)

def run_temp2_absence(scripts, bam, reference_2bit, te_bed, insert_size, threads, out, config, log):
    mccutils.log("temp2","running TEMP2 non-reference absence prediction", log=log)
    command = [
        "bash", scripts+"/TEMP2", "absence", 
            "-i", bam,
            "-r", te_bed,
            "-t", reference_2bit,
            "-f", str(insert_size),
            "-c", str(threads),
            "-o", out
    ]

    for param in config.PARAMS['absence'].keys():
        command.append(param)
        command.append(str(config.PARAMS['absence']))

    mccutils.run_command(command, log=log)

if __name__ == "__main__":                
    main()