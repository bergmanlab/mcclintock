import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.popoolationte2.popoolationte2_post as config

def main():
    mccutils.log("popoolationte2","processing PopoolationTE2 results")
    te_predictions = snakemake.input.popoolationte2_out
    te_gff = snakemake.input.te_gff
    taxonomy = snakemake.input.taxonomy
    out_dir = snakemake.params.out_dir
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    log = snakemake.params.log

    ref_tes = get_ref_tes(te_gff, taxonomy, chromosomes)
    insertions = read_insertions(te_predictions, ref_tes, chromosomes, sample_name, both_end_support_needed=config.REQUIRE_BOTH_END_SUPPORT, support_threshold=config.FREQUENCY_THRESHOLD)
    if len(insertions) >= 1:
        insertions = mccutils.make_redundant_bed(insertions, sample_name, out_dir, method="popoolationte2")
        mccutils.make_nonredundant_bed(insertions, sample_name, out_dir, method="popoolationte2")
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_popoolationte2_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_popoolationte2_nonredundant.bed"])
    
    mccutils.log("popoolationte2","PopoolationTE2 postprocessing complete")

def get_ref_tes(gff, taxon, chroms):
    ref_inserts = []
    te_family = {}
    with open(taxon,"r") as t:
        for line in t:
            split_line = line.split("\t")
            te_id = split_line[0]
            family = split_line[1]
            te_family[te_id] = family

    with open(gff, "r") as g:
        for line in g:
            if "#" not in line:
                split_line = line.split("\t")
                insert = mccutils.Insertion()
                insert.type = "reference"
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])
                insert.end = int(split_line[4])
                insert.strand = split_line[6]
                insert.family = te_family[split_line[2]]
                if insert.chromosome in chroms:
                    ref_inserts.append(insert)
    
    return ref_inserts

def read_insertions(predictions, ref_tes, chroms, sample, both_end_support_needed=True, support_threshold=0.1):
    insertions = []

    with open(predictions,"r") as tsv:
        for line in tsv:
            split_line = line.split("\t")
            insert = mccutils.Insertion()
            insert.chromosome = split_line[1]
            insert.start = int(split_line[2])
            insert.end = int(split_line[2])
            insert.strand = split_line[3]
            insert.family = split_line[4]
            insert.popoolationte2.support_type = split_line[6]
            insert.popoolationte2.frequency = float(split_line[8])

            if (insert.popoolationte2.support_type == "FR" or not both_end_support_needed) and insert.popoolationte2.frequency > support_threshold:     
                # determine if insert is a ref insert       
                for x in range(0,len(ref_tes)):
                    if ref_tes[x].start <= insert.start and insert.start <= ref_tes[x].end:
                        insert.family = ref_tes[x].family
                        insert.popoolationte2.added = ref_tes[x].popoolationte2.added
                        if not ref_tes[x].popoolationte2.added:
                            ref_tes[x].popoolationte2.added = True
                            
                        insert.type = "reference"
                        insert.start = ref_tes[x].start
                        insert.end = ref_tes[x].end
                        insert.strand = ref_tes[x].strand
                
                if insert.type == "reference":
                    insert.name = insert.family+"_reference_"+str(insert.popoolationte2.frequency)+"_"+sample+"_popoolationte2_rp_"
                else:
                    insert.type = "non-reference"
                    insert.name = insert.family+"_non-reference_"+str(insert.popoolationte2.frequency)+"_"+sample+"_popoolationte2_rp_"
                
                if not insert.popoolationte2.added:
                    insertions.append(insert)
    
    return insertions


if __name__ == "__main__":                
    main()