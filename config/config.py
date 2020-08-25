# ALL_METHODS = ["ngs_te_mapper", "relocate", "relocate2", "temp", "retroseq", "popoolationte", "popoolationte2", "te-locate", "tepid", "coverage", "trimgalore"]
# SINGLE_END_METHODS = ["ngs_te_mapper","relocate", "relocate2", "coverage", "tepid", "trimgalore"]
# MULTI_THREAD_METHODS = ["coverage", "temp", "relocate2", "ngs_te_mapper", "popoolationte", "tepid", "trimgalore"]
ALL_METHODS = ["ngs_te_mapper", "relocate", "relocate2", "temp", "retroseq", "popoolationte", "popoolationte2", "te-locate", "coverage", "trimgalore"]
SINGLE_END_METHODS = ["ngs_te_mapper","relocate", "relocate2", "coverage", "trimgalore"]
MULTI_THREAD_METHODS = ["coverage", "temp", "relocate2", "ngs_te_mapper", "popoolationte", "trimgalore"]
NO_INSTALL_METHODS = ["trimgalore"]

INPUT_DIR = "{{indir}}"
REF_NAME = "{{refname}}"
SAMPLE_NAME = "{{samplename}}"

INTERMEDIATE_PATHS = {
        'mcc_files' : INPUT_DIR,
        'reference' : INPUT_DIR+"genome_fasta/"+REF_NAME+".fasta",
        'unaugmented_reference': INPUT_DIR+"genome_fasta/"+REF_NAME+"_unaugmented.fasta",
        'masked_fasta' : INPUT_DIR+"genome_fasta/"+REF_NAME+".masked.fasta",
        'popoolationTE_ref_fasta' : INPUT_DIR+"genome_fasta/"+REF_NAME+".masked.popoolationTE.fasta",
        'telocate_ref_fasta' : INPUT_DIR+"genome_fasta/"+REF_NAME+".aug.telocate.fasta",
        'ref_2bit' : INPUT_DIR+"genome_fasta/"+REF_NAME+".aug.fasta.2bit",
        'consensus' : INPUT_DIR+"consensus_fasta/consensusTEs.fasta",
        'relocaTE_consensus' : INPUT_DIR+"consensus_fasta/formattedConsensusTEs.relocaTE.fasta",
        'locations' : INPUT_DIR+"reference_te_locations/inrefTEs.gff",
        'unaugmented_locations' : INPUT_DIR+"reference_te_locations/unaugmented_inrefTEs.gff",
        'telocate_te_gff' : INPUT_DIR+"reference_te_locations/inrefTEs_HL.gff",
        'ref_tes_bed' : INPUT_DIR+"reference_te_locations/"+REF_NAME+".ref.TEs.bed",
        'relocaTE_ref_TEs' : INPUT_DIR+"reference_te_locations/"+REF_NAME+".ref.TEs.relocaTE.gff",
        'popoolationTE_gff' : INPUT_DIR+"reference_te_locations/"+REF_NAME+".ref.TEs.popoolationTE.gff",
        'taxonomy' : INPUT_DIR+"te_taxonomy/taxonomy.tsv",
        'unaugmented_taxonomy' : INPUT_DIR+"te_taxonomy/unaugmented_taxonomy.tsv",
        'popoolationTE_taxonomy' : INPUT_DIR+"te_taxonomy/taxonomy.popoolationTE.tsv",
        'coverage_fasta' : INPUT_DIR+"coverageTEs.fasta",
        'ref_te_fasta' : INPUT_DIR+REF_NAME+".ref.TEs.fasta",
        'fq1' : INPUT_DIR+"fastq/"+SAMPLE_NAME+"_1.fq",
        'fq2' : INPUT_DIR+"fastq/"+SAMPLE_NAME+"_2.fq",
        'sam' : INPUT_DIR+"mapped_reads/"+SAMPLE_NAME+".sam",
        'bam' : INPUT_DIR+"mapped_reads/"+SAMPLE_NAME+".sorted.bam",
        'telocate_sam' : INPUT_DIR+"mapped_reads/"+SAMPLE_NAME+".telocate.sam",
        'flagstat' : INPUT_DIR+SAMPLE_NAME+".bam.flagstat",
        'median_insert_size' : INPUT_DIR+"median_insert.size",
        'repeatmasker_out' : INPUT_DIR+REF_NAME+".repeatmasker.out"
    }


RESULTS_DIR = "{{results}}"
OUT_PATHS = {
        'coverage': RESULTS_DIR+"coverage/te_depth.csv",
        'ngs_te_mapper': RESULTS_DIR+"ngs_te_mapper/"+SAMPLE_NAME+"_ngs_te_mapper_nonredundant.bed",
        'relocate': RESULTS_DIR+"relocaTE/"+SAMPLE_NAME+"_relocate_nonredundant.bed",
        'temp': RESULTS_DIR+"TEMP/"+SAMPLE_NAME+"_temp_nonredundant.bed",
        'retroseq': RESULTS_DIR+"retroseq/"+SAMPLE_NAME+"_retroseq_nonredundant.bed",
        'popoolationte': RESULTS_DIR+"popoolationTE/"+SAMPLE_NAME+"_popoolationte_nonredundant.bed",
        'popoolationte2': RESULTS_DIR+"popoolationTE2/"+SAMPLE_NAME+"_popoolationte2_nonredundant.bed",
        'te-locate': RESULTS_DIR+"te-locate/"+SAMPLE_NAME+"_telocate_nonredundant.bed",
        'trimgalore': INPUT_DIR+"fastq/"+SAMPLE_NAME+"_1.fq",
        'relocate2': RESULTS_DIR+"relocaTE2/"+SAMPLE_NAME+"_relocate2_nonredundant.bed",
        'tepid': RESULTS_DIR+"tepid/"+SAMPLE_NAME+"_tepid_nonredundant.bed"
}