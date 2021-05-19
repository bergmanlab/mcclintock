
'''
RelocaTE live:
usage:
./relocaTE.pl [-t TE_fasta_file][-g chromosome_genome_fasta][-d dir_of_fq][-e short_sample_name][-h] 

options:

**required:
-t |--te_fasta		file		fasta containing nucleotide sequences of transposable elements with 
					TSD=xxx in the desc. [no default]
-d |--fq_dir		dir		directory of paired and unpaired fastq files (paired _p1.fq & _p2.fq)
					(.fq or .fastq is acceptable)  [no default]

**recommended: 
-g |--genome_fasta	file		genome (reference) fasta file path. If not provided will only align 
					reads to TE and remove TE seq from short reads. [no default]
-e |--exper 		STR		Short sample name, will be used in the output files to create IDs for
					the insert (ex. A123) [not.given]
-o |--outdir 		STR		name for directory to contain output directories and files, will be
					created for the run (ex. 04222012_A123) [outdir_teSearch]

**optional:
-p |--parallel 		INT		Break down the single big job of relocaTE into as many smaller jobs as
					possible. The alternative (0) would be to run one after the other
					(int, 0=false or 1=true) [0] 
-q |--qsub_q 		STR		same as qsub -q option, not required [no default]
-a |--qsus_array	INT		if \'-a 1\' , create qsub PBS array jobs to run the many shell scripts
					created in the \'-a 1\' option. (see: man qsub option -t).(
					0=false or 1=true) [0] 
-w |--workingdir	dir		base working directory, needs to exist, will not be creates, full path
					required [cwd] 
-l |--len		INT		len cutoff for the TE trimmed reads to be aligned [10] 
-m |--mismatch		FRACTION	mismatch allowance for alignment to TE (ex 0.1) [0] 
-1 |--mate_1_id		STR		string to uniquely identify mate 1 paired files ex: file_p1.fq [_p1]
-2 |--mate_2_id		STR		pattern to uniquely identify mate 2 paired files ex: file_p2.fq [_p2]
-u |--unpaired_id	STR		pattern to uniquely identify unpaired files ex: file.unPaired.fq [.unPaired] 
-bm|--blat_minScore	INT		blat minScore value, used by blat in the comparison of reads to TE sequence [10]
-bt|--blat_tileSize	INT		blat tileSize value, used by blat in the comparison of reads to TE sequence  [7]
-f |--flanking_seq	INT		length of the sequence flanking the found insertion to be returned. This
					sequence is taken from the reference genome [100]
-r |--reference_ins	STR		To identify reference and shared insertions (reference and reads)
					choose option-1 or option-2. 
					option-1) (recommended) use \'-r 1\' to have RelocaTE find the location of your TE in the 
					reference.
					option-2) input the file name of a tab-delimited file containing the coordinates
					of TE insertions pre-existing in the reference sequence. [no default]
-b2 |--bowtie2	        INT             to use bowtie2 use \'-b2 1\' else for bowtie use \'-b2 0\' [0]
-h  |--help				this message


See documentation for more information. http://srobb1.github.com/RelocaTE/
'''


PARAMS = {
    '-l' : 10,
    '-m' : 0.0,
    '-bm' : 10,
    '-bt' : 7,
    '-f' : 100
}