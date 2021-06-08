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
    mccutils.log("tebreak","running tebreak")
    consensus = snakemake.input.consensus_fasta
    bam = snakemake.input.bam
    reference = snakemake.input.ref_fasta
    repeatmasker_out = snakemake.input.rm_out

    threads = snakemake.threads

    script_dir = snakemake.params.script_dir
    out_dir = snakemake.params.out_dir
    ref_name = snakemake.params.ref_name
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    status_log = snakemake.params.status_log

    try:
        ref_tes = get_ref_tes(repeatmasker_out, out_dir)
        run_tebreak(bam, consensus, reference, ref_tes, script_dir, out_dir, threads, config.PARAMS, log=log)
        mccutils.check_file_exists(snakemake.output[0])
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        with open(log,"a") as l:
            print(track, file=l)
        mccutils.log("tebreak","tebreak run failed")
        with open(status_log,"w") as l:
            l.write("FAILED\n")

        mccutils.run_command(["touch", snakemake.output[0]])
    
    mccutils.log("tebreak","tebreak run complete")

def get_ref_tes(rm_out, out):
    ref_te_file = out+"/ref_tes.txt"
    with open(rm_out, "r") as inf, open(ref_te_file, "w") as outf:
        lines = []
        for ln,line in enumerate(inf):
            if ln > 2:
                parsed_line = []
                split_line = line.split(" ")
                for val in split_line:
                    if val != "":
                        parsed_line.append(val)
                lines.append(parsed_line)
        
        for split_line in lines:
            chrom = split_line[4]
            start = split_line[5]
            end = split_line[6]
            te_type = split_line[10]
            te_fam = split_line[9]
            strand = split_line[8]

            if strand == "C":
                strand = "-"
            
            out_line = "\t".join([chrom, start, end, te_type, te_fam, strand])
            outf.write(out_line+"\n")
    
    return ref_te_file
            
def run_tebreak(bam, consensus, reference, ref_tes, script_dir, out_dir, threads, params, log=None):
    os.chdir(out_dir)
    command = ["samtools", "index", bam]
    mccutils.run_command(command, log=log)
    command = [script_dir+"/tebreak", "-b", bam, '-r', reference, '-p', str(threads), '-d', ref_tes, '-i', consensus]
    for key in params.keys():
        if params[key] != False:
            if params[key] == True:
                command.append(key)
            else:
                command.append(key)
                command.append(params[key])
    mccutils.run_command(command, log=log)

if __name__ == "__main__":                
    main()