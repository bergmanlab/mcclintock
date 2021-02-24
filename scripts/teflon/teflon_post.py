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
import scripts.output as output


def main():
    mccutils.log("teflon","TEFLoN postprocessing")

    teflon_raw = snakemake.input.teflon_out
    ref_te_bed = snakemake.input.ref_bed
    reference_fasta = snakemake.input.reference_fasta

    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    out = snakemake.output.out

    ref_tes = get_ref_tes(ref_te_bed)
    insertions = read_insertions(
        teflon_raw, 
        chromosomes, 
        sample_name, 
        ref_tes,
        min_presence=config.PARAMETERS['min_presence_reads'], 
        max_absence=config.PARAMETERS['max_absence_reads'],
        min_presence_fraction=config.PARAMETERS['min_presence_fraction'],
        require_tsd=config.PARAMETERS['require_tsd'],
        require_both_breakpoints=config.PARAMETERS['require_both_breakpoints']
    )
    if len(insertions) >= 1:
        insertions = output.make_redundant_bed(insertions, sample_name, out_dir, method="teflon")
        insertions = output.make_nonredundant_bed(insertions, sample_name, out_dir, method="teflon")
        output.write_vcf(insertions, reference_fasta, sample_name, "teflon", out_dir)
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_teflon_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_teflon_nonredundant.bed"])


def get_ref_tes(ref_te_bed):
    ref_tes = {}
    with open(ref_te_bed,"r") as bed:
        for line in bed:
            split_line = line.split("\t")
            ref_tes[split_line[3]] = [split_line[0], int(split_line[1]), int(split_line[2])]
    
    return ref_tes

def read_insertions(predictions, chroms, sample, ref_tes, min_presence=3, max_absence=None, min_presence_fraction=0.1, require_tsd=False, require_both_breakpoints=False):
    insertions = []

    with open(predictions, "r") as tsv:
        for line in tsv:
            split_line = line.split("\t")
            insert = output.Insertion(output.Teflon())

            insert.chromosome = split_line[0]

            both_ends = False
            tsd = False
            if insert.chromosome in chroms:
                if split_line[1] != "-" and split_line[2] != "-":
                    left = int(split_line[1])
                    right = int(split_line[2])
                    both_ends = True

                elif split_line[1] == "-":
                    left = int(split_line[2])
                    right = int(split_line[2])
                else:
                    left = int(split_line[1])
                    right = int(split_line[1])
                
                if left > right:
                    tsd = True
                    tmp = right
                    right = left
                    left = tmp
                
                elif left == right:
                    tsd = True
                    right += 1
                

                insert.start = left-1
                insert.end = right

                insert.family = split_line[3]
                
                insert.strand = split_line[5]
                
                # if reference prediction, uses ref TE coordinates
                if split_line[6] != "-":
                    tsd = True
                    both_ends = True
                    insert.type = "reference"
                    te_names = split_line[6].split(",")
                    te_name = ""
                    for name in te_names:
                        if name in ref_tes.keys():
                            if te_name == "":
                                te_name = name
                    
                    if te_name == "":
                        sys.exit("TEFLON ERROR: can't find:"+split_line[6]+" in reference TEs...\n")
                    insert.chromosome = ref_tes[te_name][0]
                    insert.start = ref_tes[te_name][1]
                    insert.end = ref_tes[te_name][2]

                else:
                    insert.type = "non-reference"
                
                insert.support_info.support['five_prime_supported'].value = split_line[7]
                insert.support_info.support['three_prime_supported'].value = split_line[7]
                
                insert.support_info.support['presence_reads'].value = int(split_line[9])
                insert.support_info.support['absence_reads'].value = int(split_line[10])
                insert.support_info.support['ambiguous_reads'].value = int(split_line[11])
                insert.support_info.support['frequency'].value = float(split_line[12])

                insert.name = insert.family+"|"+insert.type+"|"+str(insert.support_info.support['frequency'].value)+"|"+sample+"|teflon|rp|"

                if (
                    (insert.support_info.support['presence_reads'].value >= min_presence) 
                    and (max_absence is None or insert.support_info.support['absence_reads'].value <= max_absence)
                    and (insert.support_info.support['frequency'].value >= min_presence_fraction)
                    and ((tsd or not require_tsd) and (both_ends or not require_both_breakpoints))
                ):
                    insertions.append(insert)
    
    return insertions
            



if __name__ == "__main__":                
    main()