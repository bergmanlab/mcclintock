import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.telocate.telocate_run as config


def main():
    print("<TELOCATE> Running TE-Locate...")
    te_gff = snakemake.input.te_gff
    sam = snakemake.input.sam
    ref_fasta = snakemake.input.ref
    median_insert_size_file = snakemake.input.median_insert_size
    log = snakemake.params.log
    with open(log,"a") as l:
        l.write("TE GFF: "+te_gff+"\n")
        l.write("SAM: "+sam+"\n")
        l.write("reference fasta: "+ref_fasta+"\n")
        

    telocate = snakemake.params.run_script
    max_mem = snakemake.params.max_mem
    out_dir = snakemake.params.out_dir
    

    sam_dir = out_dir+"/sam/"
    mccutils.mkdir(sam_dir)
    te_locate_sam = sam_dir+"te-locate.sam"
    os.symlink(sam, te_locate_sam)

    os.chdir(os.path.dirname(telocate))

    median_insert_size = mccutils.get_median_insert_size(median_insert_size_file)

    distance = (median_insert_size * config.MIN_DISTANCE)

    command = ["perl", telocate, str(max_mem), sam_dir, te_gff, ref_fasta, out_dir, str(distance), str(config.MIN_SUPPORT_READS), str(config.MIN_SUPPORT_INDIVIDUALS)]

    mccutils.run_command(command, log=log)


    mccutils.run_command(["cp", out_dir+"_"+str(distance)+"_reads3_acc1.info", out_dir+"te-locate-raw.info"])
    print("<TELOCATE> TE-Locate complete")








if __name__ == "__main__":                
    main()

