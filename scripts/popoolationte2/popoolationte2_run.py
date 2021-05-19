import os
import sys
import subprocess
import importlib.util as il
import traceback
spec = il.spec_from_file_location("config", snakemake.params.config)
config = il.module_from_spec(spec)
sys.modules[spec.name] = config
spec.loader.exec_module(config)
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils


def main():
    mccutils.log("popoolationte2","running PopoolationTE2")
    ref_fasta = snakemake.input.ref_fasta
    bam = snakemake.input.bam
    taxonomy = snakemake.input.taxonomy
    jar = snakemake.params.jar
    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    log = snakemake.params.log
    status_log = snakemake.params.status_log

    prev_step_succeeded = mccutils.check_status_file(status_log)

    if prev_step_succeeded:
        try:
            mccutils.mkdir(out_dir+"/tmp")
            taxonomy = format_taxonomy(taxonomy, out_dir)
            ppileup = popoolationte2_ppileup(jar, config.PARAMS["ppileup"], bam, taxonomy, out_dir, log=log)
            ppileup = popoolationte2_subsample(jar, config.PARAMS["subsampleppileup"], ppileup, out_dir, log=log)
            signatures = popoolationte2_signatures(jar, config.PARAMS["identifySignatures"], ppileup, out_dir, log=log)
            signatures = popoolationte2_strand(jar, config.PARAMS["updateStrand"], signatures, bam, taxonomy, out_dir, log=log)
            signatures = popoolationte2_frequency(jar, ppileup, signatures, out_dir, log=log)
            te_insertions = popoolationte2_pairup(jar, config.PARAMS["pairupSignatures"], signatures, ref_fasta, taxonomy, out_dir, log=log)
            mccutils.remove(out_dir+"/tmp")
            mccutils.check_file_exists(snakemake.output[0])

            with open(status_log,"w") as l:
                l.write("COMPLETED\n")
            mccutils.log("popoolationte2","popoolationte2 run complete")
        
        except Exception as e:
            track = traceback.format_exc()
            print(track, file=sys.stderr)
            with open(log,"a") as l:
                print(track, file=l)
            mccutils.log("popoolationte2","popoolationte2 run failed")
            with open(status_log,"w") as l:
                l.write("FAILED\n")
            
            mccutils.run_command(["touch", snakemake.output[0]])
    
    else:
        mccutils.run_command(["touch", snakemake.output[0]])

def format_taxonomy(taxon, out):
    out_taxon = out+"/input.taxonomy.txt"
    with open(out_taxon,"w") as out:
        with open(taxon, "r") as tax:
            for ln, line in enumerate(tax):
                split_line = line.split("\t")
                if ln == 0:
                    out.write("id\tfamily\torder\n")
                else:
                    out.write(split_line[0]+"\t"+split_line[1]+"\tna\n")
    return out_taxon

def popoolationte2_ppileup(jar, params, bam, taxon, out, log=None):
    mccutils.log("popoolationte2","making physical pileup file", log=log)
    ppileup = out+"/output.ppileup.gz"
    command = ['java', "-Djava.io.tmpdir="+out+"/tmp", "-jar", jar, "ppileup", 
                                               "--bam", bam, 
                                               "--hier", taxon, 
                                               "--output", ppileup]
   
    
    for param in params.keys():
        command.append(param)
        command.append(str(params[param]))

    mccutils.run_command(command, log=log)

    return ppileup

def popoolationte2_subsample(jar, params, ppileup, out, log=None):
    out_ppileup = out+"/subsampled.ppileup.gz"
    if params["run"]:
        mccutils.log("popoolationte2","subsampling physical pileup file to uniform coverage", log=log)
        command = ["java", "-Djava.io.tmpdir="+out+"/tmp", "-jar", jar, "subsampleppileup",
                                       "--ppileup", ppileup,
                                       "--output", out_ppileup]

        for param in params.keys():
            if param == True and param != "run":
                command.append(param)
            elif param != False:
                command.append(param)
                command.append(str(params[param]))

        mccutils.run_command(command, log=log)
    else:
        out_ppileup = ppileup
    
    return out_ppileup


def popoolationte2_signatures(jar, params, ppileup, out, log=None):
    mccutils.log("popoolationte2","identifying signatures of TE insertions", log=log)
    signatures = out+"/output.signatures"
    command = ["java", "-Djava.io.tmpdir="+out+"/tmp", "-jar", jar, "identifySignatures",
                                               "--ppileup", ppileup,
                                               "--mode", "separate",
                                               "--output", signatures]

    for param in params.keys():
        command.append(param)
        command.append(str(params[param]))

    mccutils.run_command(command, log=log)
    
    return signatures

def popoolationte2_strand(jar, params, signatures, bam, taxon, out, log=None):
    mccutils.log("popoolationte2", "estimating strand of TEs", log=log)
    out_sig = out+"/output.stranded.signatures"
    command = ["java", "-Djava.io.tmpdir="+out+"/tmp", "-jar", jar, "updateStrand",
                                              "--signature", signatures,
                                              "--output", out_sig,
                                              "--bam", bam,
                                              "--hier", taxon]

    for param in params.keys():
        command.append(param)
        command.append(str(params[param]))

    mccutils.run_command(command, log=log)
    
    return out_sig

def popoolationte2_frequency(jar, ppileup, signatures, out, log=None):
    mccutils.log("popoolationte2","estimating frequencies for signatures of TE insertions", log=log)
    freq_signatures = out+"/output.stranded.signatures.freq"
    mccutils.run_command(["java", "-Djava.io.tmpdir="+out+"/tmp", "-jar", jar, "frequency",
                                               "--ppileup", ppileup,
                                               "--signature", signatures,
                                               "--output", freq_signatures], log=log)
    
    return freq_signatures

def popoolationte2_pairup(jar, params, signatures, ref, taxon, out, log=None):
    mccutils.log("popoolationte2","generating raw TE insertion predictions", log=log)
    te_insertions = out+"/teinsertions.txt"
    command = ["java", "-Djava.io.tmpdir="+out+"/tmp", "-jar", jar, "pairupSignatures",
                                              "--signature", signatures,
                                              "--ref-genome", ref,
                                              "--hier", taxon,
                                              "--output-detail", "medium",
                                              "--output", te_insertions]

    for param in params.keys():
        command.append(param)
        command.append(str(params[param]))

    mccutils.run_command(command, log=log)
    
    return te_insertions


if __name__ == "__main__":                
    main()