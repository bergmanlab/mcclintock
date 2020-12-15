import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import mccutils as mccutils

class Insertion:
    def __init__(self, support_info):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "None"
        self.strand = "."
        self.family = ""
        self.support_info = support_info


class Info:
    def __init__(self, tag, desc, value):
        self.tag = tag
        self.description = desc
        self.value = value

class Ngs_te_mapper:
    def __init__(self):
        self.support = {
            "supportingreads": Info("SUPPORTING_READS", "Total number of reads supporting the start and end positions", 0)
        }

class Temp:
    def __init__(self):
        self.support = {
            "class": Info(
                "CLASS", 
                """The class of the insertion. “1p1” means that the detected insertion is supported by reads at both sides. 
                “2p” means the detected insertion is supported by more than 1 read at only 1 side. 
                “Singleton” means the detected insertion is supported by only 1 read at 1 side""", 
                ""
            ),
            "variantsupport" : Info("RP_SUPPORT", "The total number of read pairs that support the detected insertion", 0),
            "frequency" : Info("FREQUENCY", "The estimated population frequency of the detected insertion", 0.0),
            "junction1" : Info("JUNCTION1", "The coordinate of the start junction. If the junction is not found, will be the arithmetic mean of the start and end coordinates", 0),
            "junction1support" : Info("JUNCTION1SUPPORT", "The number of the reads supporting the start junction. If the junction is not found, will have the value 0", 0),
            "junction2" : Info("JUNCTION2", "The coordinate of the end junction. If the junction is not found, will be the arithmetic mean of the start and end coordinates", 0),
            "junction2support" : Info("JUNCTION2SUPPORT", "The number of the reads supporting the end junction. If the junction is not found,will have the value 0", 0),
            "fiveprimesupport" : Info("FIVE_PRIME_SUPPORT", "The number of reads supporting the detected insertion at the 5’ end of the TE (not including junction spanning reads)", 0),
            "threeprimesupport" : Info("THREE_PRIME_SUPPORT", "The number of reads supporting the detected insertion at the 3’ end of the TE (not including junction spanning reads)", 0)
        }

class Telocate:
    def __init__(self):
        self.support = {
            "read_pair_support" : Info("RP_SUPPORT", "The total number of all supporting read pairs", 0)
        }

class Retroseq:
    def __init__(self):
        self.support = {
            "read_pair_support" : Info("RP_SUPPORT", "Number of correctly mapped read pairs spanning breakpoint", 0),
            "clip3" : Info("CLIP3", "Number of soft clipped reads downstream of the breakpoint", 0),
            "clip5" : Info("CLIP5", "Number of soft clipped reads upstream of the breakpoint", 0),
            "call_status" : Info(
                "CALL_STATUS", 
                """Call Status - for reference calls a flag to say if the call failed a particular filter. 
                Filters are ordered by priority in calling (higher number indicates closer to being called). 
                1 - depth too high in region, 2 - not enough reads in cluster, 3 - not enough total flanking reads, 
                4 - not enough inconsistently mapped reads, 5 - neither side passes ratio test, 6 - one side passes ratio test, 
                7 - distance too large at breakpoint, 8 - PASSED all filters""",
                0
            )
        }

class Relocate2:
    def __init__(self):
        self.support = {
            "right_junction_reads" : Info("RIGHT_JUNCTION_READS", "Number of reads covering the junction of TE insertion on right side/downstream", 0),
            "left_junction_reads" : Info("LEFT_JUNCTION_READS", "Number of reads covering the junction of TE insertion on left side/upstream", 0),
            "right_support_reads" : Info("RIGHT_SUPPORT_READS", "Number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on right side/downstream", 0),
            "left_support_reads" : Info("LEFT_SUPPORT_READS", "Number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on left side/downstream", 0)
        }

class Relocate:
    def __init__(self):
        self.support = {
            "right_flanking_reads" : Info("RIGHT_FLANKING_READS", "Number of reads that cover the right junction of the insertion site", 0),
            "left_flanking_reads" : Info("LEFT_FLANKING_READS", "Number of reads that cover the left junction of the insertion site", 0)
        }

class Popoolationte:
    def __init__(self):
        self.support = {
            "flanks_supported" : Info("FLANKS_SUPPORTED", " is the TE insertion supported by a forward (F), by a reverse (R) or by both (FR) insertions", ""),
            "frequency" : Info("FREQUENCY", "population frequency (1..fixed)", 0.0),
            "forward_insert_start" : Info("FORWARD_INSERT_START", "start of the range of the forward insertion", 0),
            "forward_insert_end" : Info("FORWARD_INSERT_END", "end of the range of the forward insertion", 0),
            "forward_insert_freq" : Info("FORWARD_INSERT_FREQ", "population frequency estimated by the forward insertion", 0.0),
            "forward_insert_cov" : Info("FORWARD_INSERT_COV", "coverage of the forward insertion", 0),
            "forward_presence_reads" : Info("FORWARD_PRESENCE_READS", "TE-presence reads of the forward insertion", 0),
            "forward_absence_reads" : Info("FORWARD_ABSENCE_READS", "TE-absence reads of the forward insertion", 0),
            "reverse_insert_start" : Info("REVERSE_INSERT_START", "start of the range of the reverse insertion", 0),
            "reverse_insert_end" : Info("REVERSE_INSERT_END", "end of the range of the reverse insertion", 0),
            "reverse_insert_freq" : Info("REVERSE_INSERT_FREQ", "population frequency estimated by the reverse insertion", 0.0),
            "reverse_insert_cov" : Info("REVERSE_INSERT_COV", "coverage of the reverse insertion", 0),
            "reverse_presence_reads" : Info("REVERSE_PRESENCE_READS", "TE-presence reads of the reverse insertion", 0),
            "reverse_absence_reads" : Info("REVERSE_ABSENCE_READS", "TE-absence reads of the reverse insertion", 0)
        }

class Popoolationte2:
    def __init__(self):
        self.support = {
            "flanks_supported" : Info("FLANKS_SUPPORTED", "support for the TE insertions; either a single forward signature (F) or a single reverse signature (R) or a matching pair of forward and reverse signatures (FR)", ""),
            "frequency" : Info("FREQUENCY", "the population frequency of the TE insertions", 0.0)
        }
        self.added = False

class Teflon:
    def __init__(self):
        self.support = {
            "five_prime_supported" : Info("FIVE_PRIME_SUPPORTED", "5' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")", ""),
            "three_prime_supported" : Info("THREE_PRIME_SUPPORTED", "3' breakpoint is supported by soft-clipped reads (if TRUE "+" else "-")", ""),
            "presence_reads" : Info("PRESENCE_READS", "read count for presence reads", 0),
            "absence_reads" : Info("ABSENCE_READS", "read count for absence reads", 0),
            "ambiguous_reads" : Info("AMBIGUOUS_READS", "read count for ambiguous reads", 0),
            "frequency" : Info("FREQUENCY", "allele frequency", 0.0)
        }

#################################################
## TODO convert to proper format for VCF creation
class Tepid:
    def __init__(self):
        self.id = -1
        self.support = 0

class Jitterbug:
    def __init__(self):
        supporting_fwd_reads = 0
        supporting_rev_reads = 0
        split_read_support = 0
        zygosity = 0.0
################################################


def make_redundant_bed(insertions, sample_name, out_dir, method="popoolationte"):
    tmp_bed = out_dir+"/tmp.bed"

    insertion_dict = {}
    out_inserts = []
    malformed_inserts = []
    properly_formed_inserts = []
    for insert in insertions:
        if insert.start <= insert.end:
            insertion_dict[ "_".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])] = insert
            properly_formed_inserts.append(insert)
        else:
            malformed_inserts.append(insert)

    insertions = properly_formed_inserts
    # write malformed predictions to separate bed file
    if len(malformed_inserts) > 0:
        malformed_bed = out_dir+"/"+sample_name+"_"+method+"_malformed.bed"
        with open(malformed_bed,"w") as out:
            for insert in malformed_inserts:
                out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
                out.write(out_line+"\n")

    with open(tmp_bed, "w") as out:
        for insert in insertions:
            out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
            out.write(out_line+"\n")
    
    sorted_bed = out_dir+"/sorted.bed"
    command = ["bedtools", "sort", "-i", tmp_bed]
    mccutils.run_command_stdout(command, sorted_bed)

    redundant_bed = out_dir+"/"+sample_name+"_"+method+"_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_'+method+'" description="'+sample_name+'_'+method+'"\n'
        outbed.write(header)
        with open(sorted_bed, "r") as inbed:
            for x, line in enumerate(inbed):

                # outputs inserts in sorted order with unique number added to name
                key = line.replace("\t","_")
                key = key.replace("\n","")
                insert = insertion_dict[key]
                insert.name += str(x+1)
                out_inserts.append(insert)

                # write to bed with unique number added to name
                split_line = line.split("\t")
                split_line[3] += str(x+1)
                line = "\t".join(split_line)
                outbed.write(line)
    
    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)

    return out_inserts

def make_nonredundant_bed(insertions, sample_name, out_dir, method="popoolationte"):
    uniq_inserts = {}

    for insert in insertions:
        key = "_".join([insert.chromosome, str(insert.start), str(insert.end), insert.type])
        if key not in uniq_inserts.keys():
            uniq_inserts[key] = insert
        else:
            ## method specific way to determine which duplicate to keep
            if method == "popoolationte":
                if (uniq_inserts[key].popoolationte.f_read_support + uniq_inserts[key].popoolationte.r_read_support) >  (insert.popoolationte.f_read_support + insert.popoolationte.r_read_support):
                    uniq_inserts[key] = insert
            
            elif method == "popoolationte2":
                if (uniq_inserts[key].popoolationte2.frequency) >  (insert.popoolationte2.frequency):
                    uniq_inserts[key] = insert
            
            elif method == "relocate":
                if (uniq_inserts[key].relocate.left_support + uniq_inserts[key].relocate.right_support) < (insert.relocate.left_support + insert.relocate.right_support):
                    uniq_inserts[key] = insert
            
            elif method == "relocate2":
                if (uniq_inserts[key].relocate2.left_support + uniq_inserts[key].relocate2.right_support) < (insert.relocate2.left_support + insert.relocate2.right_support):
                    uniq_inserts[key] = insert
            
            elif method == "retroseq":
                if insert.support_info.support['read_pair_support'].value > uniq_inserts[key].support_info.support['read_pair_support'].value:
                    uniq_inserts[key] = insert
            
            elif method == "te-locate":
                if insert.support_info.support['read_pair_support'].value > uniq_inserts[key].support_info.support['read_pair_support'].value:
                    uniq_inserts[key] = insert
            
            elif method == "temp":
                if insert.support_info.support['variantsupport'].value > uniq_inserts[key].support_info.support['variantsupport'].value:
                    uniq_inserts[key] = insert
            
            elif method == "ngs_te_mapper":
                if insert.support_info.support['supportingreads'].value > uniq_inserts[key].support_info.support['supportingreads'].value:
                    uniq_inserts[key] = insert
            
            elif method == "tepid":
                if insert.tepid.support > uniq_inserts[key].tepid.support:
                    uniq_inserts[key] = insert
            
            elif method == "teflon":
                if insert.teflon.presence_reads > uniq_inserts[key].teflon.presence_reads:
                    uniq_inserts[key] = insert
            
            elif method == "jitterbug":
                if insert.jitterbug.split_read_support > uniq_inserts[key].jitterbug.split_read_support:
                    uniq_inserts[key] = insert
    
    tmp_bed = out_dir+"/tmp.bed"
    with open(tmp_bed, "w") as outbed:
        for key in uniq_inserts.keys():
            insert = uniq_inserts[key]
            out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
            outbed.write(out_line+"\n")
    
    sorted_bed = out_dir+"/sorted.bed"
    command = ["bedtools", "sort", "-i", tmp_bed]
    mccutils.run_command_stdout(command, sorted_bed)

    nonredundant_bed = out_dir+"/"+sample_name+"_"+method+"_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_'+method+'" description="'+sample_name+'_'+method+'"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)

    out_inserts = []
    with open(nonredundant_bed, "r") as outbed:
        for line in outbed:
            if "track name=" not in line and "description=" not in line:
                split_line = line.split("\t")
                ref_type = split_line[3].split("|")[1]
                key = "_".join([split_line[0], str(int(split_line[1])+1), split_line[2], ref_type])
                out_inserts.append(uniq_inserts[key])
    mccutils.remove(tmp_bed)
    mccutils.remove(sorted_bed)

    return out_inserts

def write_vcf(inserts, genome_fasta, sample_name, method, out_dir):
    contig_lengths = {}
    with open(genome_fasta+".fai","r") as fai:
        for line in fai:
            split_line = line.split("\t")
            contig_lengths[split_line[0]] = split_line[1]
    

    contigs_with_inserts = []
    for insert in inserts:
        if insert.chromosome not in contigs_with_inserts:
            contigs_with_inserts.append(insert.chromosome)
    

    out_vcf = out_dir+"/"+sample_name+"_"+method+"_nonredundant_non-reference.vcf"
    with open(out_vcf, "w") as vcf:
        today = date.today()
        today_date = today.strftime("%Y-%m-%d")
        meta = [
            "##fileformat=VCFv4.2",
            "##fileDate="+today_date,
            "##source=McClintock",
            "##reference="+genome_fasta
        ]
        for contig in contigs_with_inserts:
            meta.append("##contig=<ID="+contig+",length="+contig_lengths[contig]+">")
        
        meta.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structure variant">')
        meta.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
        meta.append('##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation">')
        meta.append('##INFO=<ID=FAMILY,Number=1,Type=String,Description="TE family">')

        for line in meta:
            vcf.write(line+"\n")
        header = "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])
        vcf.write(header+"\n")

        ##TODO get ref sequence
        ref = "N"
        for insert in inserts:
            vals = [insert.chromosome, str(insert.start), ".", ref, "<INS:ME>", ".", "."]
            info = ["END="+str(insert.end),"SVTYPE=INS", "STRAND="+insert.strand, "FAMILY="+insert.family]

            out_line = ("\t".join(vals)) + (";".join(info))
            vcf.write(out_line+"\n")