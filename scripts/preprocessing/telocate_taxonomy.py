import os
import sys
import subprocess
import traceback
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    log = snakemake.params.log
    tmp_dir = snakemake.params.tmp_dir
    mccutils.mkdir(tmp_dir+"/telocate")

    mccutils.log("processing","making TE-locate taxonomy file", log=log)
    try:
        mccutils.run_command(["cp", snakemake.input.ref_gff, "telocate_locations.gff"])
        mccutils.run_command(["cp", snakemake.input.taxonomy, "telocate_taxonomy.tsv"])
        command = ["perl", snakemake.input.script, "telocate_locations.gff", "telocate_taxonomy.tsv", "Alias"]
        mccutils.run_command(command, log=log)
        mccutils.run_command(["cp", "telocate_locations_HL.gff", snakemake.output[0]])
        mccutils.check_file_exists(snakemake.output[0])
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...unable to produce TE-locate taxonomy file using", snakemake.input.script, file=sys.stderr)
        sys.exit(1)

    mccutils.log("processing","TE-locate taxonomy file created")

if __name__ == "__main__":                
    main()