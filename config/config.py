
ALL_METHODS = ["ngs_te_mapper", "relocate", "relocate2", "temp", "retroseq", "popoolationte", "te-locate", "coverage", "trimgalore"]
SINGLE_END_METHODS = ["ngs_te_mapper","relocate", "relocate2", "coverage", "trimgalore"]

MULTI_THREAD_METHODS = ["coverage", "temp", "relocate2", "ngs_te_mapper", "popoolationte"]

INPUT_DIR = "{{indir}}"
REF_NAME = "{{refname}}"
SAMPLE_NAME = "{{samplename}}"

INTERMEDIATE_PATHS = {
        'mcc_files' : INPUT_DIR,
        'reference' : INPUT_DIR+REF_NAME+".fasta",
        'consensus' : INPUT_DIR+"consensusTEs.fasta",
        'fq1' : INPUT_DIR+"fastq/"+SAMPLE_NAME+"_1.fq",
        'fq2' : INPUT_DIR+"fastq/"+SAMPLE_NAME+"_2.fq",
        'locations' : INPUT_DIR+"inrefTEs.gff",
        'taxonomy' : INPUT_DIR+"taxonomy.tsv",
        'coverage_fasta' : INPUT_DIR+"coverageTEs.fasta",
        'masked_fasta' : INPUT_DIR+REF_NAME+".masked.fasta",
        'formatted_ref_tes' : INPUT_DIR+REF_NAME+".ref.TEs.gff",
        'formatted_taxonomy' : INPUT_DIR+REF_NAME+".TE.taxonomy.tsv",
        'formatted_consensus' : INPUT_DIR+"formattedConsensusTEs.fasta",
        'ref_te_fasta' : INPUT_DIR+REF_NAME+".ref.TEs.fasta",
        'augmented_reference' : INPUT_DIR+REF_NAME+".aug.fasta",
        'ref_tes_bed' : INPUT_DIR+REF_NAME+".ref.TEs.bed",
        'telocate_te_gff' : INPUT_DIR+REF_NAME+".ref.TEs_HL.gff",
        'telocate_sam' : INPUT_DIR+"telocate_sam/"+SAMPLE_NAME+".telocate.sam",
        'telocate_ref_fasta' : INPUT_DIR+REF_NAME+".aug.telocate.fasta",
        'sam' : INPUT_DIR+SAMPLE_NAME+".sam",
        'bam' : INPUT_DIR+SAMPLE_NAME+".sorted.bam",
        'ref_2bit' : INPUT_DIR+REF_NAME+".aug.fasta.2bit",
        'relocaTE_consensus' : INPUT_DIR+"formattedConsensusTEs.relocaTE.fasta",
        'relocaTE_ref_TEs' : INPUT_DIR+REF_NAME+".ref.TEs.relocaTE.gff",
        'popoolationTE_ref_fasta' : INPUT_DIR+REF_NAME+".masked.popoolationTE.fasta",
        'popoolationTE_taxonomy' : INPUT_DIR+"taxonomy.popoolationTE.tsv",
        'popoolationTE_gff' : INPUT_DIR+REF_NAME+".ref.TEs.popoolationTE.gff",
        'repeatmasker_out' : INPUT_DIR+REF_NAME+".repeatmasker.out",
        'flagstat' : INPUT_DIR+SAMPLE_NAME+".bam.flagstat",
        'median_insert_size' : INPUT_DIR+"median_insert.size"
    }


RESULTS_DIR = "{{results}}"
OUT_PATHS = {
        'coverage': RESULTS_DIR+"coverage/te_depth.csv",
        'ngs_te_mapper': RESULTS_DIR+"ngs_te_mapper/"+SAMPLE_NAME+"_ngs_te_mapper_nonredundant.bed",
        'relocate': RESULTS_DIR+"relocaTE/"+SAMPLE_NAME+"_relocate_nonredundant.bed",
        'temp': RESULTS_DIR+"TEMP/"+SAMPLE_NAME+"_temp_nonredundant.bed",
        'retroseq': RESULTS_DIR+"retroseq/"+SAMPLE_NAME+"_retroseq_nonredundant.bed",
        'popoolationte': RESULTS_DIR+"popoolationTE/"+SAMPLE_NAME+"_popoolationte_nonredundant.bed",
        'te-locate': RESULTS_DIR+"te-locate/"+SAMPLE_NAME+"_telocate_nonredundant.bed",
        'trimgalore': INPUT_DIR+"fastq/"+SAMPLE_NAME+"_1.fq",
        'relocate2': RESULTS_DIR+"relocaTE2/"+SAMPLE_NAME+"_relocate2_nonredundant.bed"
}