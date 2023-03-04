import sys
import os
from Bio import SeqIO
from datetime import date
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
    def __init__(self, tag, desc, value, info_type):
        self.tag = tag
        self.description = desc
        self.value = value
        self.type = info_type

# used in redundancy filtering to determine which prediction to keep (which has most support)
class RedundancyFilter:
    def __init__(self, keys, operation):
        self.keys = keys
        self.operation = operation
    
    def get_value(self, support_info):
        value = 0
        if self.operation == "value":
            value = support_info[self.keys[0]].value
        
        elif self.operation == "sum":
            for key in self.keys:
                value += support_info[key].value
        
        return value


class Ngs_te_mapper:
    def __init__(self):
        self.support = {
            "supportingreads": Info("SUPPORTING_READS", "Total number of reads supporting the start and end positions", 0, "Integer")
        }

        self.redundancy_filter = RedundancyFilter(["supportingreads"], "value")

class Ngs_te_mapper2:
    def __init__(self):
        self.support = {
            "frequency": Info("FREQUENCY", "Estimated allele frequency", 0.0, "Float"),
            "three_prime_support" : Info("THREE_PRIME_SUPPORT", "Number of reads supporting the 3' breakpoint", 0, "Integer"),
            "five_prime_support" : Info("FIVE_PRIME_SUPPORT", "Number of reads supporting the 5' breakpoint", 0, "Integer"),
            "reference_reads": Info("REFERENCE_READS", "reads supporting the reference state at this position",0, "Integer")
        }

        self.redundancy_filter = RedundancyFilter(["three_prime_support", "five_prime_support"], "sum")

class Temp:
    def __init__(self):
        self.support = {
            "class": Info(
                "CLASS", 
                "The class of the insertion. '1p1' means that the detected insertion is supported by reads at both sides. '2p' means the detected insertion is supported by more than 1 read at only 1 side. 'Singleton' means the detected insertion is supported by only 1 read at 1 side", 
                "",
                "String"
            ),
            "variantsupport" : Info("RP_SUPPORT", "The total number of read pairs that support the detected insertion", 0, "Integer"),
            "frequency" : Info("FREQUENCY", "The estimated population frequency of the detected insertion", 0.0, "Float"),
            "junction1" : Info("JUNCTION1", "The coordinate of the start junction. If the junction is not found, will be the arithmetic mean of the start and end coordinates", 0, "Integer"),
            "junction1support" : Info("JUNCTION1SUPPORT", "The number of the reads supporting the start junction. If the junction is not found, will have the value 0", 0, "Integer"),
            "junction2" : Info("JUNCTION2", "The coordinate of the end junction. If the junction is not found, will be the arithmetic mean of the start and end coordinates", 0, "Integer"),
            "junction2support" : Info("JUNCTION2SUPPORT", "The number of the reads supporting the end junction. If the junction is not found,will have the value 0", 0, "Integer"),
            "fiveprimesupport" : Info("FIVE_PRIME_SUPPORT", "The number of reads supporting the detected insertion at the 5’ end of the TE (not including junction spanning reads)", 0, "Integer"),
            "threeprimesupport" : Info("THREE_PRIME_SUPPORT", "The number of reads supporting the detected insertion at the 3’ end of the TE (not including junction spanning reads)", 0, "Integer")
        }

        self.redundancy_filter = RedundancyFilter(["variantsupport"], "value")

class Temp2:
    def __init__(self):
        self.support = {
            "class": Info(
                "CLASS", 
                "The class of the insertion. '1p1' means that the detected insertion is supported by reads at both sides. '2p' means the detected insertion is supported by more than 1 read at only 1 side. 'singleton' means the detected insertion is supported by only 1 read at 1 side", 
                "",
                "String"
            ),
            "frequency": Info("FREQUENCY", "Frequency of the inserted transposon. It generally means what fraction of sequenced genome present this insertion", 0.0, "Float"),
            "supportreads": Info("SUPPORT_READS", "Number of reads supporting this insertion", 0.0, "Float"),
            "referencereads": Info("REF_READS", "Number of reads that do not support this insertion, AKA reference reads", 0.0, "Float"),
            "fiveprimesupport": Info("FIVE_PRIME_SUPPORT", "Number of supporting reads at 5'end of the insertion", 0.0, 'Float'),
            "threeprimesupport": Info("THREE_PRIME_SUPPORT", "Number of supporting reads at 3'end of the insertion", 0.0, "Float"),
            "reliability": Info("RELIABILITY", "Reliability of this insertion (0–100). 100 for 2p and 1p1 insertions. For singleton insertions, TEMP2 already filtered out most of the false positives but not all of them. The reliability is a percentage stand for how many singleton insertions of a specific transposon is", 0.0, "Float"),
            "fiveprimejunctionsupport": Info("FIVE_PRIME_JUNCTION_SUPPORT", "Number of supporting reads at 5'end of the insertion junction", 0.0, "Float"),
            "threeprimejunctionsupport": Info("THREE_PRIME_JUNCTION_SUPPORT", "Number of supporting reads at 3'end of the insertion junction", 0.0, "Float")
        }

        self.redundancy_filter = RedundancyFilter(["supportreads"], "value")

class Telocate:
    def __init__(self):
        self.support = {
            "read_pair_support" : Info("RP_SUPPORT", "The total number of all supporting read pairs", 0, "Integer")
        }

        self.redundancy_filter = RedundancyFilter(["read_pair_support"], "value")

class Retroseq:
    def __init__(self):
        self.support = {
            "supporting_reads" : Info("SUPPORTING_READS", "number of reads supporting the prediction", 0, "Integer"),
            "spanning_pairs" : Info("SPANNING_PAIRS", "Number of correctly mapped read pairs spanning breakpoint", 0, "Integer"),
            "clip3" : Info("CLIP3", "Number of soft clipped reads downstream of the breakpoint", 0, "Integer"),
            "clip5" : Info("CLIP5", "Number of soft clipped reads upstream of the breakpoint", 0, "Integer"),
            "frequency" : Info("FREQUENCY", "Frequency of the insertion. supporting_reads / (2 x spanning_pairs + supporting_reads)", 0.0, "Float"),
            "call_status" : Info(
                "CALL_STATUS", "Call Status - for reference calls a flag to say if the call failed a particular filter. Filters are ordered by priority in calling (higher number indicates closer to being called). 1 - depth too high in region, 2 - not enough reads in cluster, 3 - not enough total flanking reads, 4 - not enough inconsistently mapped reads, 5 - neither side passes ratio test, 6 - one side passes ratio test, 7 - distance too large at breakpoint, 8 - PASSED all filters",
                0,
                "Integer"
            )
        }

        self.redundancy_filter = RedundancyFilter(["supporting_reads"], "value")

class Relocate2:
    def __init__(self):
        self.support = {
            "right_junction_reads" : Info("RIGHT_JUNCTION_READS", "Number of reads covering the junction of TE insertion on right side/downstream", 0, "Integer"),
            "left_junction_reads" : Info("LEFT_JUNCTION_READS", "Number of reads covering the junction of TE insertion on left side/upstream", 0, "Integer"),
            "right_support_reads" : Info("RIGHT_SUPPORT_READS", "Number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on right side/downstream", 0, "Integer"),
            "left_support_reads" : Info("LEFT_SUPPORT_READS", "Number of reads not covering the junction of TE insertion, but supporting TE insertion by paired-end reads on left side/downstream", 0, "Integer")
        }

        self.redundancy_filter = RedundancyFilter(["left_support_reads", "right_support_reads"], "sum")

class Relocate:
    def __init__(self):
        self.support = {
            "right_flanking_reads" : Info("RIGHT_FLANKING_READS", "Number of reads that cover the right junction of the insertion site", 0, "Integer"),
            "left_flanking_reads" : Info("LEFT_FLANKING_READS", "Number of reads that cover the left junction of the insertion site", 0, "Integer")
        }

        self.redundancy_filter = RedundancyFilter(["left_flanking_reads", "right_flanking_reads"], "sum")

class Popoolationte:
    def __init__(self):
        self.support = {
            "flanks_supported" : Info("FLANKS_SUPPORTED", "is the TE insertion supported by a forward (F), by a reverse (R) or by both (FR) insertions", "", "String"),
            "frequency" : Info("FREQUENCY", "population frequency (1..fixed)", 0.0, "Float"),
            "forward_insert_start" : Info("FORWARD_INSERT_START", "start of the range of the forward insertion", 0, "Integer"),
            "forward_insert_end" : Info("FORWARD_INSERT_END", "end of the range of the forward insertion", 0, "Integer"),
            "forward_insert_freq" : Info("FORWARD_INSERT_FREQ", "population frequency estimated by the forward insertion", 0.0, "Float"),
            "forward_insert_cov" : Info("FORWARD_INSERT_COV", "coverage of the forward insertion", 0, "Integer"),
            "forward_presence_reads" : Info("FORWARD_PRESENCE_READS", "TE-presence reads of the forward insertion", 0, "Integer"),
            "forward_absence_reads" : Info("FORWARD_ABSENCE_READS", "TE-absence reads of the forward insertion", 0, "Integer"),
            "reverse_insert_start" : Info("REVERSE_INSERT_START", "start of the range of the reverse insertion", 0, "Integer"),
            "reverse_insert_end" : Info("REVERSE_INSERT_END", "end of the range of the reverse insertion", 0, "Integer"),
            "reverse_insert_freq" : Info("REVERSE_INSERT_FREQ", "population frequency estimated by the reverse insertion", 0.0, "Float"),
            "reverse_insert_cov" : Info("REVERSE_INSERT_COV", "coverage of the reverse insertion", 0, "Integer"),
            "reverse_presence_reads" : Info("REVERSE_PRESENCE_READS", "TE-presence reads of the reverse insertion", 0, "Integer"),
            "reverse_absence_reads" : Info("REVERSE_ABSENCE_READS", "TE-absence reads of the reverse insertion", 0, "Integer")
        }

        self.redundancy_filter = RedundancyFilter(["forward_presence_reads", "reverse_presence_reads"], "sum")

class Popoolationte2:
    def __init__(self):
        self.support = {
            "flanks_supported" : Info("FLANKS_SUPPORTED", "support for the TE insertions; either a single forward signature (F) or a single reverse signature (R) or a matching pair of forward and reverse signatures (FR)", "", "String"),
            "frequency" : Info("FREQUENCY", "the population frequency of the TE insertions", 0.0, "Float")
        }
        self.added = False

        self.redundancy_filter = RedundancyFilter(["frequency"], "value")

class Teflon:
    def __init__(self):
        self.support = {
            "five_prime_supported" : Info("FIVE_PRIME_SUPPORTED", "5' breakpoint is supported by soft-clipped reads (if TRUE '+' else '-')", "", "String"),
            "three_prime_supported" : Info("THREE_PRIME_SUPPORTED", "3' breakpoint is supported by soft-clipped reads (if TRUE '+' else '-')", "", "String"),
            "presence_reads" : Info("PRESENCE_READS", "read count for presence reads", 0, "Integer"),
            "absence_reads" : Info("ABSENCE_READS", "read count for absence reads", 0, "Integer"),
            "ambiguous_reads" : Info("AMBIGUOUS_READS", "read count for ambiguous reads", 0, "Integer"),
            "frequency" : Info("FREQUENCY", "allele frequency", 0.0, "Float")
        }

        self.redundancy_filter = RedundancyFilter(["presence_reads"], "value")

class Tebreak:
    def __init__(self):
        self.support = {
            "five_p_elt_match" : Info("FIVE_P_ELT_MATCH", "Fraction of bases matched to reference for inserted sequence on insertion seqment of 5' supporting contig", 0.0, "Float"),
            "three_p_elt_match" : Info("THREE_P_ELT_MATCH", "Fraction of bases matched to reference for inserted sequence on insertion seqment of 3' supporting contig", 0.0, "Float"),
            "five_p_genome_match" : Info("FIVE_P_GENOME_MATCH", "Fraction of bases matched to reference genome on genomic segment of 5' supporting contig", 0.0, "Float"),
            "three_p_genome_match" : Info("THREE_P_GENOME_MATCH", "Fraction of bases matched to reference genome on genomic segment of 3' supporting contig", 0.0, "Float"),
            "split_reads_5prime" : Info("SPLIT_READS_5PRIME", "Number of split reads supporting 5' end of the insertion", 0, "Integer"),
            "split_reads_3prime" : Info("SPLIT_READS_3PRIME", "Number of split reads supporting 3' end of the insertion", 0, "Integer"),
            "remapped_discordant" : Info("REMAPPED_DISCORDANT", "Number of discordant read ends re-mappable to insertion reference sequence", 0, "Integer"),
            "remap_disc_fraction" : Info("REMAP_DISC_FRACTION", "The proportion of remapped discordant reads mapping to the reference insertion sequence", 0.0, "Float"),
            "remapped_splitreads" : Info("REMAPPED_SPLITREADS", "Number of split reads re-mappable to insertion reference sequence", 0, "Integer"),
            "remap_split_fraction" : Info("REMAP_SPLIT_FRACTION", "The proportion of remapped split reads mapping to the reference insertion sequence", 0.0, "Float")
        }

        self.redundancy_filter = RedundancyFilter(["split_reads_5prime", "split_reads_3prime"], "sum")

class Jitterbug:
    def __init__(self):
        self.support = {
            "supporting_fwd_reads" : Info("SUPPORTING_FWD_READS", "Reads supporting the forward junction of the insertion", 0, "Integer"),
            "supporting_rev_reads" : Info("SUPPORTING_REV_READS", "Reads supporting the reverse junction of the insertion", 0, "Integer"),
            "softclipped_support" : Info("SOFTCLIPPED_SUPPORT", "Number of softclipped reads that support the breakpoint", 0, "Integer"),
            "zygosity" : Info("ZYGOSITY", "The ratio of clipped and properly mapped reads at the insertion site", 0.0, "Float")
        }

        self.redundancy_filter = RedundancyFilter(["supporting_fwd_reads", "supporting_rev_reads"], "sum")

class Tepid:
    def __init__(self):
        self.id = -1
        self.support = {
            "supporting_reads" : Info("SUPPORTING_READS", "Number of reads supporting the insertion", 0, "Integer")
        }
        self.redundancy_filter = RedundancyFilter(["supporting_reads"], "value")


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
            if (insert.support_info.redundancy_filter.get_value(insert.support_info.support) > uniq_inserts[key].support_info.redundancy_filter.get_value(uniq_inserts[key].support_info.support)):
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

def write_vcf(inserts, genome_fasta, sample_name, method, out_dir, vcf_options):
    contigs = {}

    for record in SeqIO.parse(genome_fasta, "fasta"):
        seq_name = str(record.id)
        seq = str(record.seq)
        contigs[seq_name] = seq
    

    contigs_with_inserts = []
    for insert in inserts:
        if insert.chromosome not in contigs_with_inserts:
            contigs_with_inserts.append(insert.chromosome)

    # create site-only vcf   
    if "siteonly" in vcf_options:
        out_vcf = out_dir+"/"+sample_name+"_"+method+"_nonredundant_non-reference_siteonly.vcf"
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
                meta.append("##contig=<ID="+contig+",length="+str(len(contigs[contig]))+">")
            
            meta.append('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structure variant">')
            meta.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
            meta.append('##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation">')
            meta.append('##INFO=<ID=FAMILY,Number=1,Type=String,Description="TE family">')
            for key,value in inserts[0].support_info.support.items():
                meta.append('##INFO=<ID='+value.tag+',Number=1,Type='+value.type+',Description="'+value.description+'">')

            for line in meta:
                vcf.write(line+"\n")
            header = "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])
            vcf.write(header+"\n")


            for insert in inserts:
                if insert.type == "non-reference":
                    ref = contigs[insert.chromosome][insert.start-1]
                    te_id = insert.family+"_"+insert.name.split("|")[-1]
                    vals = [insert.chromosome, str(insert.start), te_id, ref.upper(), "<INS:ME>", ".", "PASS"]
                    info = ["END="+str(insert.end),"SVTYPE=INS", "STRAND="+insert.strand, "FAMILY="+insert.family]
                    for key,value in insert.support_info.support.items():
                        info.append(value.tag+"="+str(value.value))

                    out_line = ("\t".join(vals)) + "\t" + (";".join(info))
                    vcf.write(out_line+"\n")

    # create mergable vcf with sample column
    if "sample" in vcf_options:
        out_vcf_sample = out_dir+"/"+sample_name+"_"+method+"_nonredundant_non-reference_sample.vcf"
        with open(out_vcf_sample, "w") as vcf:
            today = date.today()
            today_date = today.strftime("%Y-%m-%d")
            meta = [
                "##fileformat=VCFv4.2",
                "##fileDate="+today_date,
                "##source=McClintock",
                "##reference="+genome_fasta
            ]
            for contig in contigs_with_inserts:
                meta.append("##contig=<ID="+contig+",length="+str(len(contigs[contig]))+">")
            
            meta.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
            meta.append('##INFO=<ID=STRAND,Number=1,Type=String,Description="Strand orientation">')
            meta.append('##INFO=<ID=FAMILY,Number=1,Type=String,Description="TE family">')
            meta.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
            meta.append('##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position of the structure variant">')
            for key,value in inserts[0].support_info.support.items():
                meta.append('##FORMAT=<ID='+value.tag+',Number=1,Type='+value.type+',Description="'+value.description+'">')

            for line in meta:
                vcf.write(line+"\n")
            header = "\t".join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sample_name])
            vcf.write(header+"\n")


            for insert in inserts:
                if insert.type == "non-reference":
                    ref = contigs[insert.chromosome][insert.start-1]
                    te_id = insert.family+"_"+insert.name.split("|")[-1]
                    vals = [insert.chromosome, str(insert.start), te_id, ref.upper(), "<INS:ME>", ".", "PASS"]
                    info = ["SVTYPE=INS", "STRAND="+insert.strand, "FAMILY="+insert.family]
                    format_col = ["GT","END"]
                    sample_col = ["1/1",str(insert.end)]
                    for key,value in insert.support_info.support.items():
                        format_col.append(value.tag)
                        sample_col.append(str(value.value))

                    out_line = ("\t".join(vals)) + "\t" + (";".join(info)) + "\t" + (":".join(format_col)) + "\t" + (":".join(sample_col))
                    vcf.write(out_line+"\n")