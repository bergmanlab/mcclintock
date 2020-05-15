import os
import sys
import subprocess
import traceback
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    log = snakemake.params.log
    print("<PROCESSING> making TE-locate taxonomy file...log:"+log)
    try:
        command = ["perl", snakemake.input.script, snakemake.input.ref_gff, snakemake.input.taxonomy, "Alias"]
        mccutils.run_command(command, log=log)
        mccutils.check_file_exists(snakemake.output[0])
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...unable to produce TE-locate taxonomy file using", snakemake.input.script, file=sys.stderr)
        sys.exit(1)

        

if __name__ == "__main__":                
    main()