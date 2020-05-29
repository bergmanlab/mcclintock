import os
import sys
import subprocess
from Bio import SeqIO
import traceback
try:
    sys.path.append(snakemake.config['args']['mcc_path'])
    import scripts.mccutils as mccutils
    import scripts.fix_fasta as fix_fasta
except Exception as e:
    track = traceback.format_exc()
    print(track, file=sys.stderr)
    print("ERROR...unable to locate required external scripts at: "+snakemake.config['args']['mcc_path']+"/scripts/", file=sys.stderr)
    sys.exit(1)


def main():
    reference = snakemake.input.ref
    consensus = snakemake.input.consensus
    te_gff = snakemake.params.te_gff
    taxonomy = snakemake.params.taxonomy
    mcc_out = snakemake.params.mcc_out
    run_id = snakemake.params.run_id
    log = snakemake.params.log
    chromosomes = snakemake.params.chromosomes.split(",")
    augment = snakemake.params.augment
    processors = snakemake.threads
    out_te_gff = snakemake.output.te_gff
    out_taxonomy = snakemake.output.taxonomy

    if not os.path.exists(mcc_out+"/tmp"):
        mccutils.mkdir(mcc_out+"/tmp")

    if te_gff == "None" and taxonomy == "None":
        print("<PROCESSING> no reference TEs provided... finding reference TEs with RepeatMasker...")
        te_gff = repeat_mask(reference, consensus, chromosomes, processors, run_id, log, mcc_out)
        taxonomy, te_gff = make_te_taxonomy_map(te_gff, consensus, run_id, mcc_out)
    
    taxonomy = mccutils.replace_special_chars_taxonomy(taxonomy, mcc_out+"/tmp/"+str(run_id)+"taxonomy.tsv")
    te_gff = format_ref_te_gff(te_gff, run_id, mcc_out)

    if augment != "None":
        augment = mccutils.replace_special_chars_fasta(augment, mcc_out+"/tmp/"+str(run_id)+"augment.tmp")
        taxonomy = augment_taxonomy(taxonomy, augment, mcc_out+"/tmp/"+str(run_id)+"taxonomy2.tsv")
        te_gff = augment_gff(te_gff, augment, mcc_out+"/tmp/"+str(run_id)+"reference_tes.gff")
    
    mccutils.run_command(["mv", te_gff, out_te_gff])
    mccutils.run_command(["mv", taxonomy, out_taxonomy])





# creates reference TE gff using repeatmasker and consensus TE fasta
def repeat_mask(reference, te_fasta, chromosomes, procs, run_id, log, out):
    try:
        outdir = out+"/tmp/repeatmasker_"+run_id
        mccutils.mkdir(outdir)
        os.chdir(outdir)
        command = ["RepeatMasker","-pa", str(procs), "-lib", te_fasta, "-dir", outdir, "-s", "-gff", "-nolow", "-no_is", reference]
        mccutils.run_command(command, log=log)
        os.chdir(out)

        # RepeatMasker appears to override the custom database names during the ProcessRepeats
        # this step changes them back, more rules may be needed for other reference genomes
        ref_name = os.path.basename(reference)
        repeatmasker_gff = outdir+"/"+ref_name+".out.gff"
        formatted_ref_tes = out+"/tmp/"+run_id+"tmpreferenceTEs.gff"
        with open(repeatmasker_gff,"r") as rmgff:
            with open(formatted_ref_tes,"w") as outgff:
                for line in rmgff:
                    if "#" not in line:
                        line = line.replace("McClintock-int","McClintock")
                        line = line.replace("POGON1", "pogo")
                        split_line = line.split("\t")
                        feats = split_line[8]
                        if split_line[0] in chromosomes:
                            te = feats.split(" ")[1]
                            te = te.replace('"','').split(":")[1]
                            feats = ";".join(["ID="+te, "Name="+te, "Alias="+te])
                            split_line[2] = te
                            split_line[8] = feats
                            line = "\t".join(split_line)

                            outgff.write(line+'\n')


        masked_fasta = outdir+"/"+ref_name+".masked"
        fasta_lines = fix_fasta.fix_fasta_lines(masked_fasta, 80)

        mccutils.check_file_exists(formatted_ref_tes)

    except Exception as e:
        track = traceback.format_exc()
        print(track, file=sys.stderr)
        print("ERROR...Failed to run repeatmasker on: ", reference, "with lib:", te_fasta, "check file formatting...exiting...", file=sys.stderr)
        sys.exit(1)

    return formatted_ref_tes


def make_te_taxonomy_map(ref_tes, consensus, run_id, out):
    consensus_families = []
    with open(consensus,"r") as fa:
        for line in fa:
            if ">" in line:
                line = line.replace("\n","")
                family = line.replace(">","")
                if family not in consensus_families:
                    consensus_families.append(family)

    te_family_counts = {}
    family_tsv = out+"/tmp/"+run_id+"tmpfam.tsv"
    gff_lines = []
    with open(ref_tes, "r") as ingff:
        with open(family_tsv,"w") as outtsv:
            for line in ingff:
                if "#" not in line:
                    split_line = line.split("\t")
                    family = split_line[2]
                    if family in te_family_counts.keys():
                        te_family_counts[family] += 1
                    else:
                        te_family_counts[family] = 1
                    
                    te_id = family+"_"+str(te_family_counts[family])
                    split_line[2] = te_id
                    split_line[8] = ";".join(["ID="+te_id, "Name="+te_id, "Alias="+te_id])
                    line = "\t".join(split_line)
                    gff_lines.append(line+"\n")
                    if family in consensus_families:
                        outtsv.write(te_id+"\t"+family+"\n")
    
    with open(ref_tes, "w") as outgff:
        for line in gff_lines:
            outgff.write(line)

    return family_tsv, ref_tes

# standardizes format of reference te gff
def format_ref_te_gff(ref_tes, run_id, out):
    format_ref_tes = out+"/tmp/"+run_id+"tmpreferenceTEs1.gff"
    with open(ref_tes,"r") as ingff:
        with open(format_ref_tes,"w") as outgff:
            for line in ingff:
                if "#" not in line:
                    split_line = line.split("\t")
                    features = split_line[8].replace("\n","")
                    split_feats = features.split(";")
                    te_id = "MISSING"
                    for feat in split_feats:
                        if "ID=" in feat:
                            te_id = feat.split("=")[1]
                    
                    te_id = mccutils.replace_special_chars(te_id)
                    split_line[2] = te_id
                    features = ";".join(["ID="+te_id, "Name="+te_id, "Alias="+te_id])
                    line = "\t".join(split_line[0:8])
                    line = line+"\t"+features+"\n"
                    outgff.write(line)
    
    return format_ref_tes

def augment_taxonomy(taxonomy, fasta, out):
    lines = []
    with open(taxonomy, "r") as taxon:
        for line in taxon:
            lines.append(line)
    
    with open(fasta, "r") as fa:
        for line in fa:
            if ">" in line:
                element = line.replace(">", "")
                element = element.replace("\n","")
                line = "mcc"+element+"\t"+element+"\n"
                lines.append(line)
    
    with open(out, "w") as o:
        for line in lines:
            o.write(line)
    
    return out

def augment_gff(gff, fasta, out):
    lines = []
    with open(gff, "r") as g:
        for line in g:
            lines.append(line)
    
    fasta_records = SeqIO.parse(fasta,"fasta")
    for record in fasta_records:
        length = len(str(record.seq))
        te = str(record.id)
        features = ";".join(["ID=mcc"+te,"Name=mcc"+te,"Alias=mcc"+te])
        line = "\t".join([te, "reannotate", "mcc"+te, "1", str(length), ".", "+", ".", features])
        line = line+"\n"
        lines.append(line)

    with open(out, "w") as o:
        for line in lines:
            o.write(line)
    
    return out

if __name__ == "__main__":
    main()