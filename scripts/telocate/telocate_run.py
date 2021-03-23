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
    te_gff = snakemake.input.te_gff
    sam = snakemake.input.sam
    ref_fasta = snakemake.input.ref
    median_insert_size_file = snakemake.input.median_insert_size
    log = snakemake.params.log
    status_log = snakemake.params.status_log

    mccutils.log("te-locate","running TE-Locate", log=log)
    with open(log,"a") as l:
        l.write("TE GFF: "+te_gff+"\n")
        l.write("SAM: "+sam+"\n")
        l.write("reference fasta: "+ref_fasta+"\n")
        

    telocate = snakemake.params.run_script
    out_dir = snakemake.params.out_dir

    try:
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

        mccutils.check_file_exists(out_dir+"_"+str(distance)+"_reads3_acc1.info")
        mccutils.run_command(["cp", out_dir+"_"+str(distance)+"_reads3_acc1.info", out_dir+"te-locate-raw.info"])
    
        mccutils.log("te-locate", "TE-Locate complete")
        with open(status_log,"w") as l:
            l.write("COMPLETED\n")

    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        with open(log,"a") as l:
            print(track, file=l)
        mccutils.log("telocate","TE-locate run failed")
        with open(status_log,"w") as l:
            l.write("FAILED\n")
        
        mccutils.run_command(["touch", snakemake.output[0]])

    



if __name__ == "__main__":                
    main()

