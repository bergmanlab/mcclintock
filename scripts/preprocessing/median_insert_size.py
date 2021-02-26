import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import statistics


def main():
    fq2 = snakemake.params.fq2

    if fq2 != "None":
        mccutils.log("processing","calculating median insert size of reads")
        median = mccutils.calc_median_insert_size(snakemake.input[0])

        if median > 0:
            with open(snakemake.output[0],"w") as out:
                out.write("median_insert_size="+str(median)+"\n")
                
            mccutils.log("processing","median insert size of reads calculated")
    
    else:
        with open(snakemake.output[0],"w") as out:
            out.write("median_insert_size=0\n")
        

if __name__ == "__main__":                
    main()