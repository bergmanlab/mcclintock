ALL_METHODS = ["ngs_te_mapper", "ngs_te_mapper2", "relocate", "relocate2", "temp", "retroseq", "popoolationte", "popoolationte2", "te-locate", "teflon", "coverage", "trimgalore","map_reads"]
SINGLE_END_METHODS = ["ngs_te_mapper", "ngs_te_mapper2", "relocate", "coverage", "trimgalore"] # relocate2 removed due to bugs
MULTI_THREAD_METHODS = ["coverage", "temp", "relocate2", "ngs_te_mapper", "ngs_te_mapper2", "popoolationte", "teflon", "trimgalore"]
NO_INSTALL_METHODS = ["trimgalore", "map_reads"]

INPUT_DIR = "{{indir}}"
REF_DIR = "{{refdir}}"
SAM_DIR = "{{samdir}}"
REF_NAME = "{{refname}}"
SAMPLE_NAME = "{{samplename}}"

INTERMEDIATE_PATHS = {
        'mcc_files' : REF_DIR,
        'reference' : REF_DIR+"genome_fasta/"+REF_NAME+".fasta",
        'unaugmented_reference': REF_DIR+"genome_fasta/"+REF_NAME+"_unaugmented.fasta",
        'masked_fasta' : SAM_DIR+"intermediate/genome_fasta/"+REF_NAME+".masked.fasta",
        'popoolationTE_ref_fasta' : SAM_DIR+"intermediate/genome_fasta/"+REF_NAME+".masked.popoolationTE.fasta",
        'telocate_ref_fasta' : SAM_DIR+"intermediate/genome_fasta/"+REF_NAME+".aug.telocate.fasta",
        'ref_2bit' : SAM_DIR+"intermediate/genome_fasta/"+REF_NAME+".aug.fasta.2bit",
        'consensus' : REF_DIR+"consensus_fasta/consensusTEs.fasta",
        'relocaTE_consensus' : SAM_DIR+"intermediate/consensus_fasta/formattedConsensusTEs.relocaTE.fasta",
        'locations' : REF_DIR+"reference_te_locations/inrefTEs.gff",
        'unaugmented_locations' : REF_DIR+"reference_te_locations/unaugmented_inrefTEs.gff",
        'telocate_te_gff' : SAM_DIR+"intermediate/reference_te_locations/inrefTEs_HL.gff",
        'ref_tes_bed' : SAM_DIR+"intermediate/reference_te_locations/"+REF_NAME+".ref.TEs.bed",
        'relocaTE_ref_TEs' : SAM_DIR+"intermediate/reference_te_locations/"+REF_NAME+".ref.TEs.relocaTE.gff",
        'popoolationTE_gff' : SAM_DIR+"intermediate/reference_te_locations/"+REF_NAME+".ref.TEs.popoolationTE.gff",
        'taxonomy' : REF_DIR+"te_taxonomy/taxonomy.tsv",
        'unaugmented_taxonomy' : REF_DIR+"te_taxonomy/unaugmented_taxonomy.tsv",
        'popoolationTE_taxonomy' : SAM_DIR+"intermediate/te_taxonomy/taxonomy.popoolationTE.tsv",
        'coverage_fasta' : SAM_DIR+"intermediate/coverageTEs.fasta",
        'ref_te_fasta' : SAM_DIR+"intermediate/"+REF_NAME+".ref.TEs.fasta",
        'fq1' : SAM_DIR+"intermediate/fastq/"+SAMPLE_NAME+"_1.fq",
        'fq2' : SAM_DIR+"intermediate/fastq/"+SAMPLE_NAME+"_2.fq",
        'sam' : SAM_DIR+"intermediate/mapped_reads/"+SAMPLE_NAME+".sam",
        'bam' : SAM_DIR+"intermediate/mapped_reads/"+SAMPLE_NAME+".sorted.bam",
        'telocate_sam' : SAM_DIR+"intermediate/mapped_reads/"+SAMPLE_NAME+".telocate.sam",
        'flagstat' : SAM_DIR+"intermediate/mapped_reads/"+SAMPLE_NAME+".bam.flagstat",
        'median_insert_size' : SAM_DIR+"intermediate/mapped_reads/median_insert.size",
        'repeatmasker_out' : SAM_DIR+"intermediate/"+REF_NAME+".repeatmasker.out"
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
        'trimgalore': SAM_DIR+"intermediate/fastq/"+SAMPLE_NAME+"_1.fq",
        'map_reads': SAM_DIR+"intermediate/mapped_reads/median_insert.size",
        'relocate2': RESULTS_DIR+"relocaTE2/"+SAMPLE_NAME+"_relocate2_nonredundant.bed",
        'tepid': RESULTS_DIR+"tepid/"+SAMPLE_NAME+"_tepid_nonredundant.bed",
        'teflon': RESULTS_DIR+"teflon/"+SAMPLE_NAME+"_teflon_nonredundant.bed",
        'jitterbug': RESULTS_DIR+"jitterbug/"+SAMPLE_NAME+"_jitterbug_nonredundant.bed"
}

ESSENTIAL_PATHS = {
        'coverage': [
                RESULTS_DIR+"coverage/te_depth.csv",
                RESULTS_DIR+"coverage/plots/",
                RESULTS_DIR+"coverage/te-depth-files/"
        ],

        'ngs_te_mapper': [
                RESULTS_DIR+"ngs_te_mapper/"+SAMPLE_NAME+"_ngs_te_mapper_nonredundant.bed",
                RESULTS_DIR+"ngs_te_mapper/"+SAMPLE_NAME+"_ngs_te_mapper_redundant.bed",
                RESULTS_DIR+"ngs_te_mapper/"+SAMPLE_NAME+"_ngs_te_mapper_malformed.bed",
                RESULTS_DIR+"ngs_te_mapper/"+SAMPLE_NAME+"_ngs_te_mapper_nonredundant_non-reference.vcf",
                RESULTS_DIR+"ngs_te_mapper/unfiltered/bed_tsd/"+SAMPLE_NAME+"_insertions.bed"
        ],

        'relocate': [
                RESULTS_DIR+"relocaTE/"+SAMPLE_NAME+"_relocate_nonredundant.bed",
                RESULTS_DIR+"relocaTE/"+SAMPLE_NAME+"_relocate_redundant.bed",
                RESULTS_DIR+"relocaTE/"+SAMPLE_NAME+"_relocate_malformed.bed",
                RESULTS_DIR+"relocaTE/"+SAMPLE_NAME+"_relocate_nonredundant_non-reference.vcf",
                RESULTS_DIR+"relocaTE/unfiltered/combined.gff"
        ],

        'temp': [
                RESULTS_DIR+"TEMP/"+SAMPLE_NAME+"_temp_nonredundant.bed",
                RESULTS_DIR+"TEMP/"+SAMPLE_NAME+"_temp_redundant.bed",
                RESULTS_DIR+"TEMP/"+SAMPLE_NAME+"_temp_malformed.bed",
                RESULTS_DIR+"TEMP/"+SAMPLE_NAME+".absent.bed",
                RESULTS_DIR+"TEMP/"+SAMPLE_NAME+"_temp_nonredundant_non-reference.vcf",
                RESULTS_DIR+"TEMP/unfiltered/"+SAMPLE_NAME+".insertion.refined.bp.summary",
                RESULTS_DIR+"TEMP/unfiltered/"+SAMPLE_NAME+".absence.refined.bp.summary"
        ],

        'retroseq': [
                RESULTS_DIR+"retroseq/"+SAMPLE_NAME+"_retroseq_nonredundant.bed",
                RESULTS_DIR+"retroseq/"+SAMPLE_NAME+"_retroseq_redundant.bed",
                RESULTS_DIR+"retroseq/"+SAMPLE_NAME+"_retroseq_malformed.bed",
                RESULTS_DIR+"retroseq/"+SAMPLE_NAME+"_retroseq_nonredundant_non-reference.vcf",
                RESULTS_DIR+"retroseq/unfiltered/"+SAMPLE_NAME+".call"
        ],

        'popoolationte': [
                RESULTS_DIR+"popoolationTE/"+SAMPLE_NAME+"_popoolationte_nonredundant.bed",
                RESULTS_DIR+"popoolationTE/"+SAMPLE_NAME+"_popoolationte_redundant.bed",
                RESULTS_DIR+"popoolationTE/"+SAMPLE_NAME+"_popoolationte_malformed.bed",
                RESULTS_DIR+"popoolationTE/"+SAMPLE_NAME+"_popoolationte_nonredundant_non-reference.vcf",
                RESULTS_DIR+"popoolationTE/unfiltered/te-poly-filtered.txt"
        ],

        'popoolationte2': [
                RESULTS_DIR+"popoolationTE2/"+SAMPLE_NAME+"_popoolationte2_nonredundant.bed",
                RESULTS_DIR+"popoolationTE2/"+SAMPLE_NAME+"_popoolationte2_redundant.bed",
                RESULTS_DIR+"popoolationTE2/"+SAMPLE_NAME+"_popoolationte2_malformed.bed",
                RESULTS_DIR+"popoolationTE2/"+SAMPLE_NAME+"_popoolationte2_nonredundant_non-reference.vcf",
                RESULTS_DIR+"popoolationTE2/unfiltered/teinsertions.txt"
        ],

        'te-locate': [
                RESULTS_DIR+"te-locate/"+SAMPLE_NAME+"_telocate_nonredundant.bed",
                RESULTS_DIR+"te-locate/"+SAMPLE_NAME+"_telocate_redundant.bed",
                RESULTS_DIR+"te-locate/"+SAMPLE_NAME+"_telocate_malformed.bed",
                RESULTS_DIR+"te-locate/"+SAMPLE_NAME+"_telocate_nonredundant_non-reference.vcf",
                RESULTS_DIR+"te-locate/unfiltered/te-locate-raw.info"
        ],

        'trimgalore': [],

        'map_reads': [],

        'relocate2': [
                RESULTS_DIR+"relocaTE2/"+SAMPLE_NAME+"_relocate2_nonredundant.bed",
                RESULTS_DIR+"relocaTE2/"+SAMPLE_NAME+"_relocate2_redundant.bed",
                RESULTS_DIR+"relocaTE2/"+SAMPLE_NAME+"_relocate2_malformed.bed",
                RESULTS_DIR+"relocaTE2/"+SAMPLE_NAME+"_relocate2_nonredundant_non-reference.vcf",
                RESULTS_DIR+"relocaTE2/unfiltered/repeat/results/ALL.all_nonref_insert.gff",
                RESULTS_DIR+"relocaTE2/unfiltered/repeat/results/ALL.all_ref_insert.gff"
        ],

        'tepid': [
                RESULTS_DIR+"tepid/"+SAMPLE_NAME+"_tepid_nonredundant.bed",
                RESULTS_DIR+"tepid/"+SAMPLE_NAME+"_tepid_redundant.bed",
                RESULTS_DIR+"tepid/"+SAMPLE_NAME+"_tepid_malformed.bed",
                RESULTS_DIR+"tepid/unfiltered/insertions_"+REF_NAME+".bed",
                RESULTS_DIR+"tepid/unfiltered/deletions_"+REF_NAME+".bed",
                RESULTS_DIR+"tepid/unfiltered/insertion_reads_"+REF_NAME+".txt",
                RESULTS_DIR+"tepid/unfiltered/deletion_reads_"+REF_NAME+".txt"
        ],

        'teflon': [
                RESULTS_DIR+"teflon/"+SAMPLE_NAME+"_teflon_nonredundant.bed",
                RESULTS_DIR+"teflon/"+SAMPLE_NAME+"_teflon_redundant.bed",
                RESULTS_DIR+"teflon/"+SAMPLE_NAME+"_teflon_malformed.bed",
                RESULTS_DIR+"teflon/"+SAMPLE_NAME+"_teflon_nonredundant_non-reference.vcf",
                RESULTS_DIR+"teflon/unfiltered/genotypes/sample.genotypes.txt",
                RESULTS_DIR+"teflon/unfiltered/reference_te.bed"
        ],

        'jitterbug': [
                RESULTS_DIR+"jitterbug/"+SAMPLE_NAME+"_jitterbug_nonredundant.bed",
                RESULTS_DIR+"jitterbug/"+SAMPLE_NAME+"_jitterbug_redundant.bed",
                RESULTS_DIR+"jitterbug/"+SAMPLE_NAME+"_jitterbug_malformed.bed",
                RESULTS_DIR+"jitterbug/unfiltered/"+SAMPLE_NAME+".TE_insertions_paired_clusters.filtered.gff3"
        ]
}