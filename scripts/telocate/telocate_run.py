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
    te_gff = snakemake.input.te_gff
    sam = snakemake.input.sam
    ref_fasta = snakemake.input.ref
    median_insert_size_file = snakemake.input.median_insert_size
    log = snakemake.params.log

    mccutils.log("te-locate","running TE-Locate", log=log)
    with open(log,"a") as l:
        l.write("TE GFF: "+te_gff+"\n")
        l.write("SAM: "+sam+"\n")
        l.write("reference fasta: "+ref_fasta+"\n")
        

    telocate = snakemake.params.run_script
    out_dir = snakemake.params.out_dir

    # ensures intermediate files from previous runs are removed
    for f in os.listdir(out_dir):
        mccutils.remove(out_dir+"/"+f)

    sam_dir = out_dir+"/sam/"
    mccutils.mkdir(sam_dir)
    te_locate_sam = sam_dir+"te-locate.sam"
    if os.path.exists(te_locate_sam):
        os.remove(te_locate_sam)
    os.symlink(sam, te_locate_sam)

    os.chdir(os.path.dirname(telocate))

    median_insert_size = mccutils.get_median_insert_size(median_insert_size_file)

    distance = (median_insert_size * config.MIN_DISTANCE)

    command = ["perl", telocate, str(config.MAX_MEM), sam_dir, te_gff, ref_fasta, out_dir, str(distance), str(config.MIN_SUPPORT_READS), str(config.MIN_SUPPORT_INDIVIDUALS)]

    mccutils.run_command(command, log=log)


    mccutils.run_command(["cp", out_dir+"_"+str(distance)+"_reads3_acc1.info", out_dir+"te-locate-raw.info"])
    mccutils.log("te-locate", "TE-Locate complete")



if __name__ == "__main__":                
    main()

