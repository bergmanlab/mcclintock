ALL_METHODS = ["ngs_te_mapper", "ngs_te_mapper2", "relocate", "relocate2", "temp", "temp2", "retroseq", "popoolationte", "popoolationte2", "te-locate", "teflon", "coverage", "trimgalore","map_reads", "tebreak"]
SINGLE_END_METHODS = ["ngs_te_mapper", "ngs_te_mapper2", "relocate", "coverage", "trimgalore", "map_reads", "tebreak"]
MULTI_THREAD_METHODS = ["coverage", "temp", "temp2", "relocate2", "ngs_te_mapper", "ngs_te_mapper2", "popoolationte", "teflon", "trimgalore", "tebreak"]
NO_INSTALL_METHODS = ["trimgalore", "map_reads", "coverage", "relocate2"] # no source code to install for these methods, just envs

INPUT_DIR = "{{indir}}"
REF_DIR = "{{refdir}}"
SAM_DIR = "{{samdir}}"
REF_NAME = "{{refname}}"
SAMPLE_NAME = "{{samplename}}"


CONFIGS = {
        "coverage" : ["coverage/coverage.py"],
        "trimgalore" : ["trimgalore/trimgalore.py"],
        "ngs_te_mapper": ["ngs_te_mapper/ngs_te_mapper_run.py","ngs_te_mapper/ngs_te_mapper_post.py"],
        "ngs_te_mapper2":  ["ngs_te_mapper2/ngs_te_mapper2_run.py","ngs_te_mapper2/ngs_te_mapper2_post.py"],
        "popoolationte": ["popoolationte/popoolationte_run.py", "popoolationte/popoolationte_post.py"],
        "popoolationte2": ["popoolationte2/popoolationte2_run.py", "popoolationte2/popoolationte2_post.py"],
        "relocate": ["relocate/relocate_run.py", "relocate/relocate_post.py"],
        "relocate2": ["relocate2/relocate2_run.py", "relocate2/relocate2_post.py"],
        "retroseq": ["retroseq/retroseq_run.py", "retroseq/retroseq_post.py"],
        "teflon": ["teflon/teflon_run.py", "teflon/teflon_post.py"],
        "te-locate": ["telocate/telocate_run.py", "telocate/telocate_post.py"],
        "temp": ["temp/temp_run.py", "temp/temp_post.py"],
        "temp2": ["temp2/temp2_run.py", "temp2/temp2_post.py"],
        "tebreak": ["tebreak/tebreak_run.py", "tebreak/tebreak_post.py"]
}

# rules to re-run if specific config files change
CONFIG_RULES = {
        "coverage" : ["coverage"],
        "trimgalore" : ["setup_reads"],
        "ngs_te_mapper": ["ngs_te_mapper_run","ngs_te_mapper_post"],
        "ngs_te_mapper2":  ["ngs_te_mapper2_run","ngs_te_mapper2_post"],
        "popoolationte": ["popoolationTE_run", "popoolationTE_post"],
        "popoolationte2": ["popoolationTE2_run", "popoolationTE2_post"],
        "relocate": ["relocaTE_run", "relocaTE_post"],
        "relocate2": ["relocaTE2_run", "relocaTE2_post"],
        "retroseq": ["retroseq_run", "retroseq_post"],
        "teflon": ["teflon_run", "teflon_post"],
        "te-locate": ["telocate_run", "telocate_post"],
        "temp": ["run_temp", "process_temp"],
        "temp2": ["run_temp2", "process_temp2"],
        "tebreak": ["tebreak_run", "tebreak_post"]
}

LOG_DIR = "{{logdir}}"
STATUS_FILES ={
        "coverage" : LOG_DIR+"status/coverage.status",
        "ngs_te_mapper": LOG_DIR+"status/ngs_te_mapper.status",
        "ngs_te_mapper2":  LOG_DIR+"status/ngs_te_mapper2.status",
        "popoolationte": LOG_DIR+"status/popoolationte.status",
        "popoolationte2": LOG_DIR+"status/popoolationte2.status",
        "relocate": LOG_DIR+"status/relocate.status",
        "relocate2": LOG_DIR+"status/relocate2.status",
        "retroseq": LOG_DIR+"status/retroseq.status",
        "teflon": LOG_DIR+"status/teflon.status",
        "te-locate": LOG_DIR+"status/telocate.status",
        "temp": LOG_DIR+"status/temp.status",
        "temp2": LOG_DIR+"status/temp2.status",
        "tebreak": LOG_DIR+"status/tebreak.status"
}


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
        'ngs_te_mapper2_ref_TEs' : SAM_DIR+"intermediate/reference_te_locations/"+REF_NAME+".ref.TEs.ngs_te_mapper2.gff",
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
        'duplicate_maked_bam' : SAM_DIR+"intermediate/mapped_reads/"+SAMPLE_NAME+".sorted.duplicate_marked.bam",
        'telocate_sam' : SAM_DIR+"intermediate/mapped_reads/"+SAMPLE_NAME+".telocate.sam",
        'flagstat' : SAM_DIR+"intermediate/mapped_reads/"+SAMPLE_NAME+".bam.flagstat",
        'median_insert_size' : SAM_DIR+"intermediate/mapped_reads/median_insert.size",
        'repeatmasker_out' : SAM_DIR+"intermediate/"+REF_NAME+".repeatmasker.out"
    }


RESULTS_DIR = "{{results}}"


OUT_DIRS = {
        'coverage': RESULTS_DIR+"coverage/",
        'ngs_te_mapper': RESULTS_DIR+"ngs_te_mapper/",
        'ngs_te_mapper2': RESULTS_DIR+"ngs_te_mapper2/",
        'relocate': RESULTS_DIR+"relocate/",
        'temp': RESULTS_DIR+"temp/",
        'temp2': RESULTS_DIR+"temp2/",
        'retroseq': RESULTS_DIR+"retroseq/",
        'popoolationte': RESULTS_DIR+"popoolationte/",
        'popoolationte2': RESULTS_DIR+"popoolationte2/",
        'te-locate': RESULTS_DIR+"te-locate/",
        'trimgalore': SAM_DIR+"intermediate/fastq/",
        'map_reads': SAM_DIR+"intermediate/mapped_reads/",
        'relocate2': RESULTS_DIR+"relocate2/",
        'tepid': RESULTS_DIR+"tepid/",
        'teflon': RESULTS_DIR+"teflon/",
        'jitterbug': RESULTS_DIR+"jitterbug/",
        'tebreak': RESULTS_DIR+"tebreak/"
}

METHOD_DIR = "{{method}}"

OUT_PATHS = {
        'coverage': METHOD_DIR+"te_depth.csv",
        'ngs_te_mapper': METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper_nonredundant.bed",
        'ngs_te_mapper2': METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper2_nonredundant.bed",
        'relocate': METHOD_DIR+SAMPLE_NAME+"_relocate_nonredundant.bed",
        'temp': METHOD_DIR+SAMPLE_NAME+"_temp_nonredundant.bed",
        'temp2': METHOD_DIR+SAMPLE_NAME+"_temp2_nonredundant.bed",
        'retroseq': METHOD_DIR+SAMPLE_NAME+"_retroseq_nonredundant.bed",
        'popoolationte': METHOD_DIR+SAMPLE_NAME+"_popoolationte_nonredundant.bed",
        'popoolationte2': METHOD_DIR+SAMPLE_NAME+"_popoolationte2_nonredundant.bed",
        'te-locate': METHOD_DIR+SAMPLE_NAME+"_telocate_nonredundant.bed",
        'trimgalore': METHOD_DIR+SAMPLE_NAME+"_1.fq",
        'map_reads': METHOD_DIR+"median_insert.size",
        'relocate2': METHOD_DIR+SAMPLE_NAME+"_relocate2_nonredundant.bed",
        'tepid': METHOD_DIR+SAMPLE_NAME+"_tepid_nonredundant.bed",
        'teflon': METHOD_DIR+SAMPLE_NAME+"_teflon_nonredundant.bed",
        'jitterbug': METHOD_DIR+SAMPLE_NAME+"_jitterbug_nonredundant.bed",
        'tebreak': METHOD_DIR+SAMPLE_NAME+"_tebreak_nonredundant.bed"
}

ESSENTIAL_PATHS = {
        'coverage': [
                METHOD_DIR+"te_depth.csv",
                METHOD_DIR+"plots/",
                METHOD_DIR+"te-depth-files/"
        ],

        'ngs_te_mapper': [
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/bed_tsd/"+SAMPLE_NAME+"_insertions.bed"
        ],

        'ngs_te_mapper2': [
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper2_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper2_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper2_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_ngs_te_mapper2_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".nonref.bed",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".ref.bed"
        ],

        'relocate': [
                METHOD_DIR+SAMPLE_NAME+"_relocate_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_relocate_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_relocate_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_relocate_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/combined.gff"
        ],

        'temp': [
                METHOD_DIR+SAMPLE_NAME+"_temp_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_temp_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_temp_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+".absent.bed",
                METHOD_DIR+SAMPLE_NAME+"_temp_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".insertion.refined.bp.summary",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".absence.refined.bp.summary"
        ],

        'temp2': [
                METHOD_DIR+SAMPLE_NAME+"_temp2_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_temp2_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_temp2_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+".absent.bed",
                METHOD_DIR+SAMPLE_NAME+"_temp2_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".insertion.bed",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".absence.refined.bp.summary"
        ],

        'retroseq': [
                METHOD_DIR+SAMPLE_NAME+"_retroseq_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_retroseq_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_retroseq_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_retroseq_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".call"
        ],

        'popoolationte': [
                METHOD_DIR+SAMPLE_NAME+"_popoolationte_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_popoolationte_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_popoolationte_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_popoolationte_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/te-poly-filtered.txt"
        ],

        'popoolationte2': [
                METHOD_DIR+SAMPLE_NAME+"_popoolationte2_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_popoolationte2_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_popoolationte2_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_popoolationte2_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/teinsertions.txt"
        ],

        'te-locate': [
                METHOD_DIR+SAMPLE_NAME+"_telocate_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_telocate_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_telocate_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_telocate_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/te-locate-raw.info"
        ],

        'trimgalore': [],

        'map_reads': [],

        'relocate2': [
                METHOD_DIR+SAMPLE_NAME+"_relocate2_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_relocate2_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_relocate2_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_relocate2_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/repeat/results/ALL.all_nonref_insert.gff",
                METHOD_DIR+"unfiltered/repeat/results/ALL.all_ref_insert.gff"
        ],

        'tepid': [
                METHOD_DIR+SAMPLE_NAME+"_tepid_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_tepid_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_tepid_malformed.bed",
                METHOD_DIR+"unfiltered/insertions_"+REF_NAME+".bed",
                METHOD_DIR+"unfiltered/deletions_"+REF_NAME+".bed",
                METHOD_DIR+"unfiltered/insertion_reads_"+REF_NAME+".txt",
                METHOD_DIR+"unfiltered/deletion_reads_"+REF_NAME+".txt"
        ],

        'teflon': [
                METHOD_DIR+SAMPLE_NAME+"_teflon_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_teflon_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_teflon_malformed.bed",
                METHOD_DIR+SAMPLE_NAME+"_teflon_nonredundant_non-reference.vcf",
                METHOD_DIR+"unfiltered/genotypes/sample.genotypes.txt",
                METHOD_DIR+"unfiltered/reference_te.bed"
        ],

        'jitterbug': [
                METHOD_DIR+SAMPLE_NAME+"_jitterbug_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_jitterbug_redundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_jitterbug_malformed.bed",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".TE_insertions_paired_clusters.filtered.gff3"
        ],

        'tebreak': [
                METHOD_DIR+SAMPLE_NAME+"_tebreak_nonredundant.bed",
                METHOD_DIR+SAMPLE_NAME+"_tebreak_redundant.bed",
                METHOD_DIR+"unfiltered/"+SAMPLE_NAME+".tebreak.table.txt"
        ]
}