import os
import sys
import subprocess
import traceback
try:
    sys.path.append(snakemake.config['args']['mcc_path'])
    import scripts.mccutils as mccutils
except Exception as e:
    track = traceback.format_exc()
    print(track, file=sys.stderr)
    print("ERROR...unable to locate required external scripts at: "+snakemake.config['args']['mcc_path']+"/scripts/", file=sys.stderr)
    sys.exit(1)


def main():
    mccutils.log("processing","mapping reads to reference", log=snakemake.log[0])

    try:
        command = ["bwa","mem"]
        if eval(snakemake.config['args']['save_comments']):
            command.append("-C")
        
        command += ["-t", str(snakemake.threads), "-M", "-R", "@RG\\tID:"+snakemake.params.sample+"\\tSM:"+snakemake.params.sample, snakemake.input.ref, snakemake.input.fq1]

        if snakemake.config['in']['fq2'] != "None":
            command.append(snakemake.input.fq2)
        
        mccutils.run_command_stdout(command, snakemake.output[0], log=snakemake.log[0])
        
        mccutils.check_file_exists(snakemake.output[0])
    
    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        if snakemake.config['in']['fq2'] == "None":
            print("ERROR...unable to map reads (bwa mem) using reference fasta:",snakemake.input.ref,"and reads:", snakemake.input.fq1, file=sys.stderr)
        else:
            print("ERROR...unable to map reads (bwa mem) using reference fasta:",snakemake.input.ref,"and reads:", snakemake.input.fq1, snakemake.input.fq2, file=sys.stderr)
        sys.exit(1)
        
    mccutils.log("processing","read mapping complete")


if __name__ == "__main__":                
    main()