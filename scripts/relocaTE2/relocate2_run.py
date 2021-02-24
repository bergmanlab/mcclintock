import os
import sys
import subprocess
import importlib
spec = importlib.util.spec_from_file_location("config", snakemake.params.config)
config = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = config
spec.loader.exec_module(config)
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils



def main():
    sample_name = snakemake.params.sample_name
    threads = snakemake.threads
    out_dir = snakemake.params.out_dir
    median_insert_size_file = snakemake.input.median_insert_size
    log = snakemake.params.log

    # ensures intermediate files from previous runs are removed
    for f in os.listdir(out_dir):
        mccutils.remove(out_dir+"/"+f)

    is_paired = True
    if snakemake.params.raw_fq2 == "None":
        is_paired = False

    input_dir = snakemake.params.out_dir+"/input/"
    mccutils.remove(input_dir)
    mccutils.mkdir(input_dir)
    fq_dir = snakemake.params.out_dir+"/input/fastq/"
    mccutils.mkdir(fq_dir)

    reference = input_dir+"reference.fasta"
    te_seqs = input_dir+"consensus.fasta"
    rm_out = input_dir+"repeatmasker.out"

    os.symlink(snakemake.input.reference, reference)
    os.symlink(snakemake.input.te_seqs, te_seqs)
    os.symlink(snakemake.input.rm_out, rm_out)

    if is_paired:
        fq1 = fq_dir+sample_name+"_1.fq"
        fq2 = fq_dir+sample_name+"_2.fq"
        os.symlink(snakemake.input.fq1, fq1)
        os.symlink(snakemake.input.fq2, fq2)
    else:
        fq1 = fq_dir+sample_name+".unPaired.fq"
        os.symlink(snakemake.input.fq1, fq1)



    median_insert_size = get_median_insert_size(median_insert_size_file)
    output = subprocess.Popen(["which", "relocaTE2.py"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    script = output.stdout.read()
    script = script.decode()
    script = script.replace("\n","")




    mccutils.log("relocate2","running RelocaTE2", log=log)
    command =  [
        "python2", script, 
                        "-t", te_seqs,
                        "-g", reference,
                        "-r", rm_out,
                        "-o", out_dir,
                        "-s", str(median_insert_size),
                        "--run",
                        "-v", "4",
                        "-c", str(threads),
                        "--aligner", config.RELOCATE2["aligner"],
                        "--len_cut_match", str(config.RELOCATE2["len_cut_match"]),
                        "--len_cut_trim", str(config.RELOCATE2["len_cut_trim"]),
                        "--mismatch", str(config.RELOCATE2["mismatch"]),
                        "--mismatch_junction", str(config.RELOCATE2["mismatch_junction"]),
                        "-d", fq_dir
    ]

    if is_paired:
        command += ["-1", "_1", "-2", "_2"]
    
    else:
        command += ["-u", ".unPaired"]

    mccutils.run_command(command, log=log)

    mccutils.log("relocate2","RelocaTE2 run complete")


def get_median_insert_size(infile):
    median_insert_size = 0
    with open(infile,"r") as inf:
        for line in inf:
            insert = line.split("=")[1]
            insert = insert.replace("\n","")
            median_insert_size = int(float(insert))
    
    return median_insert_size

if __name__ == "__main__":                
    main()