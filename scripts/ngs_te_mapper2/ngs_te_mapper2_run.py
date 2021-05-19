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
    consensus_fasta = snakemake.input.consensus_fasta
    reference_fasta = snakemake.input.reference_fasta
    fastq1 = snakemake.input.fastq1
    fastq2 = snakemake.input.fastq2
    locations = snakemake.input.locations

    log = snakemake.params.log
    with open(log,"a") as l:
        l.write("consensus fasta: "+consensus_fasta+"\n")
        l.write("reference fasta: "+reference_fasta+"\n")
        l.write("fastq1: "+fastq1+"\n")
        l.write("fastq2: "+fastq2+"\n")


    threads = snakemake.threads
    sample_name = snakemake.params.sample_name
    script_dir = snakemake.params.script_dir
    out_dir = snakemake.params.out_dir
    status_log = snakemake.params.status_log
    out_bed_nonref = snakemake.output[0]
    out_bed_ref = snakemake.output[1]


    try:
        # ensures intermediate files from previous runs are removed
        for f in os.listdir(out_dir):
            mccutils.remove(out_dir+"/"+f)

        is_paired = True
        if snakemake.params.raw_fq2 == "None":
            is_paired = False
        
        command = [
            'python', script_dir+"/ngs_te_mapper2.py", 
                "-r", reference_fasta, 
                "-l", consensus_fasta, 
                "-t", str(threads), 
                "-o", out_dir, 
                "--keep_files",
                "-p", sample_name,
                "-a", locations,
        ]

        for key in config.PARAMS.keys():
            command.append(key)
            command.append(str(config.PARAMS[key]))

        command.append("-f")
        if is_paired:
            command.append(fastq1+","+fastq2)
        else:
            command.append(fastq1)
        
        mccutils.log("ngs_te_mapper2","running ngs_te_mapper2", log=log)
        mccutils.run_command(command, log=log)
        mccutils.check_file_exists(out_bed_ref)
        mccutils.check_file_exists(out_bed_nonref)
        with open(status_log,"w") as l:
            l.write("COMPLETED\n")

        mccutils.log("ngs_te_mapper2","ngs_te_mapper2 run complete", log=log)
        mccutils.log("ngs_te_mapper2","ngs_te_mapper2 run complete")
    
    except Exception as e:
        mccutils.log("ngs_te_mapper2","ngs_te_mapper2 run failed", log=log)
        mccutils.log("ngs_te_mapper2","ngs_te_mapper2 run failed")

        track = traceback.format_exc()
        print(track, file=sys.stderr)
        with open(log,"a") as l:
            print(track, file=l)

        with open(status_log,"w") as l:
            l.write("FAILED\n")

        mccutils.run_command(["touch", out_bed_ref])
        mccutils.run_command(["touch", out_bed_nonref])


    



if __name__ == "__main__":                
    main()