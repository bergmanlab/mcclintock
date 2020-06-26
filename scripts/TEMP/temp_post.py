import os
import sys
import subprocess
sys.path.append(snakemake.config['args']['mcc_path'])
import scripts.mccutils as mccutils
import config.TEMP.temp_post as config


class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.direction = "None"
        self.classification = "None"
        self.support = -1
        self.frequency = -1
        self.junction1 = -1
        self.junction1Support = -1
        self.junction2 = -1
        self.junction2Support = -1
        self.fivePrimeSupport = -1
        self.threePrimeSupport = -1
        self.type = "None"

def main():
    mccutils.log("temp","running TEMP post processing")
    insert_summary = snakemake.input.insert_summary
    absence_summary = snakemake.input.absence_summary
    te_gff = snakemake.input.te_gff
    log = snakemake.params.log
    sample_name = snakemake.params.sample_name
    chromosomes = snakemake.params.chromosomes.split(",")
    out_dir = snakemake.params.out_dir

    insertions = read_insertion_summary(insert_summary, sample_name)
    absence_bed = make_absence_bed(absence_summary, sample_name, out_dir)
    non_absent_ref_insertions = get_non_absent_ref_tes(te_gff, absence_bed, sample_name, out_dir, log)
    insertions += non_absent_ref_insertions
    insertions = filter_insertions(insertions, chromosomes)
    insertions = sort_insertions(insertions)
    if len(insertions) > 0:
        make_raw_bed(insertions, sample_name, out_dir)
        make_redundant_bed(insertions, sample_name, out_dir, log, acceptable_classes=config.ACCEPTABLE_INSERTION_SUPPORT_CLASSES, frequency_theshold=config.FREQUENCY_THRESHOLD)
        make_nonredundant_bed(insertions, sample_name, out_dir, log, acceptable_classes=config.ACCEPTABLE_INSERTION_SUPPORT_CLASSES, frequency_theshold=config.FREQUENCY_THRESHOLD)
    else:
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp_raw.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp_redundant.bed"])
        mccutils.run_command(["touch", out_dir+"/"+sample_name+"_temp_nonredundant.bed"])
    mccutils.log("temp","TEMP postprocessing complete")






def read_insertion_summary(infile, sample):
    insertions = []
    with open(infile,"r") as inf:
        for x,line in enumerate(inf):
            if x > 0:
                insert = Insertion()
                split_line = line.split("\t")
                if len(split_line) == 14:
                    insert.chromosome = split_line[0]
                    insert.start = int(split_line[1])-1
                    insert.end = int(split_line[2])
                    insert.name = split_line[3]+"_non-reference_"+split_line[7]+"_"+sample+"_temp_"

                    if  "antisense" in split_line[4]:
                        insert.direction = "-"
                    else:
                        insert.direction = "+"
                        
                    insert.classification = split_line[5]
                    insert.support = float(split_line[6])
                    insert.frequency = float(split_line[7])
                    insert.junction1 = int(split_line[8])
                    insert.junction1Support = int(split_line[9])
                    insert.junction2 = int(split_line[10])
                    insert.junction2Support = int(split_line[11])
                    insert.fivePrimeSupport = float(split_line[12])
                    insert.threePrimeSupport = float(split_line[13].replace("\n",""))
                    insert.type = "non-reference"

                    if insert.end >= insert.start and insert.end > 0 and insert.start > -1:

                        # if split read, use junction positions as start and end
                        if insert.junction1Support > 0 and insert.junction2Support > 0:
                            insert.start = insert.junction1 - 1
                            insert.end = insert.junction2
                            insert.name = insert.name+"sr_"

                        # read pair
                        else:
                            insert.name = insert.name+"rp_" 

                        insertions.append(insert)
                    else:
                        print("<TEMP POST> Omitting malformed line from insertion summary results:", line)
                else:
                    print("<TEMP POST> Omitting malformed line from insertion summary results:", line)
    
    return insertions


def make_absence_bed(summary_file, sample, out):
    out_bed = out+"/"+sample+".absent.bed"
    lines = []
    with open(summary_file, "r") as inf:
        for x,line in enumerate(inf):
            if x > 0:
                split_line = line.split("\t")
                new_line = "\t".join([split_line[0], split_line[1], split_line[2]])
                new_line += "\n"
                lines.append(new_line)
    
    if len(lines) < 1:
        lines.append("empty\t0\t1\n")
    
    with open(out_bed,"w") as bed:
        for line in lines:
            bed.write(line)
    
    return out_bed


def get_non_absent_ref_tes(te_gff, absence_bed, sample, out, log):
    insertions = []
    tmp_gff = out+"/tmp.ref_nonabs.gff"
    command = ["bedtools", "subtract", "-A", "-a", te_gff, "-b", absence_bed]
    mccutils.run_command_stdout(command, tmp_gff, log=log)

    with open(tmp_gff,"r") as gff:
        for line in gff:
            if "#" not in line:
                line = line.replace(";","\t")
                split_line = line.split("\t")
                insert = Insertion()
                insert.chromosome = split_line[0]
                insert.start = int(split_line[3])-1
                insert.end = int(split_line[4])
                insert.support = "!"
                insert.name = split_line[9].split("=")[1]+"_reference_"+sample+"_temp_nonab_"
                insert.direction = split_line[6]
                insert.classification = "!"
                insert.junction1Support = "!"
                insert.junction2Support = "!"
                insert.junction1 = '!'
                insert.junction2 = "!"
                insert.frequency = "!"
                insert.type = "reference"
                
                insertions.append(insert)
    
    mccutils.remove(tmp_gff)

    return insertions
    

def filter_insertions(insertions, chromosomes):
    out = []
    for insert in insertions:
        if insert.chromosome in chromosomes:
            out.append(insert)
    
    return out

def sort_insertions(insertions):
    # sorts insertions by chromosome and start position

    chromosomes = []
    insert_dict = {}
    sorted_insertions = []
    for insert in insertions:
        # adds chromosome to non-redundant list of chromosomes
        if insert.chromosome not in chromosomes:
            chromosomes.append(insert.chromosome)

        # separates insertions by chromosome, adds inserts to list based on start position
        if insert.chromosome not in insert_dict.keys():
            insert_dict[insert.chromosome] = []
            insert_dict[insert.chromosome].append(insert)
        else:
            inserted = False
            if insert.start <= insert_dict[insert.chromosome][0].start:
                insert_dict[insert.chromosome] = [insert] + insert_dict[insert.chromosome]
                inserted = True
            
            elif insert.start >= insert_dict[insert.chromosome][-1].start:
                insert_dict[insert.chromosome].append(insert)
                inserted = True

            else:
                for i in range(0,len(insert_dict[insert.chromosome])-1):
                    if insert_dict[insert.chromosome][i].start <= insert.start and insert.start <= insert_dict[insert.chromosome][i+1].start:
                        insert_dict[insert.chromosome].insert(i+1,insert)
                        inserted = True
                        break
            
            if inserted == False:
                sys.exit("<ERROR> <TEMP POST> An error occured in sorting of the insertions... exiting...\n")

    insert_num = 0
    chromosomes.sort()
    for chrom in chromosomes:
        for insert in insert_dict[chrom]:
            insert.name += str(insert_num)
            insert_num += 1
            sorted_insertions.append(insert)

    return sorted_insertions
    


def make_raw_bed(insertions, sample, out):
    raw_bed = out+"/"+sample+"_temp_raw.bed"
    
    with open(raw_bed,"w") as bed:
        header = 'track name="%s_TEMP" description="%s_TEMP"' % (sample, sample)
        bed.write(header+"\n")

        for insert in insertions:
            line = "\t".join([insert.chromosome, str(insert.start), str(insert.end), insert.name, "0", insert.direction])
            line += "\n"
            bed.write(line)



def make_redundant_bed(insertions, sample, out, log, acceptable_classes=["1p1"], frequency_theshold=0.1):
    redundant_bed = out+"/"+sample+"_temp_redundant.bed"

    with open(redundant_bed,"w") as bed:
        header = 'track name="%s_TEMP" description="%s_TEMP"' % (sample, sample)
        bed.write(header+"\n")
        for insert in insertions:
            if insert.type == "reference" or (insert.classification in acceptable_classes and insert.frequency > frequency_theshold):
                line = "\t".join([insert.chromosome, str(insert.start), str(insert.end), insert.name, "0", insert.direction])
                line += "\n"
                bed.write(line)


        
def make_nonredundant_bed(insertions, sample, out, log, acceptable_classes=["1p1"], frequency_theshold=0.1):
    unsorted_nonredundant_bed = out+"/"+sample+"_temp_unsorted_nonredundant.bed"

    collaped_insertions = {}

    # collapsing all insterts that share the same chromosome and end position (and pass thresholds)
    for insert in insertions:
        if insert.type == "reference" or (insert.classification in acceptable_classes and insert.frequency > frequency_theshold):
            if insert.type == "reference":
                # reference TEs are only considered 'redundant' if they share the same start and end
                key = insert.chromosome+"_"+str(insert.start)+"_"+str(insert.end)
            else:
                key = insert.chromosome+"_"+str(insert.end)

            if key not in collaped_insertions.keys():
                collaped_insertions[key] = []

            collaped_insertions[key].append(insert)
    
    with open(unsorted_nonredundant_bed, "w") as bed:
        for key in collaped_insertions.keys():
            highest_supported_insert = None
            for x, insert in enumerate(collaped_insertions[key]):
                if x < 1:
                    highest_supported_insert = insert
                else:
                    if highest_supported_insert.support != "!" and insert.support > highest_supported_insert.support:
                        highest_supported_insert = insert
            
            line = "\t".join([highest_supported_insert.chromosome, str(highest_supported_insert.start), str(highest_supported_insert.end), highest_supported_insert.name, "0", highest_supported_insert.direction])
            bed.write(line+"\n")

    tmp_bed = out+"/"+sample+"_temp_nonredundant.bed.tmp"

    command = ["bedtools", "sort", "-i", unsorted_nonredundant_bed]
    mccutils.run_command_stdout(command, tmp_bed, log=log)

    nonredundant_bed = out+"/"+sample+"_temp_nonredundant.bed"
    with open(nonredundant_bed,"w") as outbed:
        with open(tmp_bed, "r") as inbed:
            header = 'track name="%s_TEMP" description="%s_TEMP"' % (sample, sample)
            outbed.write(header+"\n")
            for line in inbed:
                outbed.write(line)

    mccutils.remove(tmp_bed)
    mccutils.remove(unsorted_nonredundant_bed)






if __name__ == "__main__":                
    main()