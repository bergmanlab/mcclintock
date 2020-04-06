
# TEMP_Insertion.sh
'''
Options:
        -i     Input file in bam format with full path. Please sort and index the file before calling this program. 
               Sorting and indexing can be done by 'samtools sort' and 'samtools index'
        -s     Directory where all the scripts are
        -o     Path to output directory. Default is current directory
        -r     Transposon sequence database in fasta format with full path
        -t     Annotated TEs in BED6 format with full path. Detected insertions that overlap with annoated TEs will be filtered. 
        -u     TE families annotations. If supplied detected insertions overlap with annotated TE of the same family will be filtered. Only use with -t.
        -m     Number of mismatch allowed when mapping to TE concensus sequences. Default is 3
        -x     The minimum score difference between the best hit and the second best hit for considering a read as uniquely mapped. For BWA mem. 
        -f     An integer specifying the length of the fragments (inserts) of the library. Default is 500
        -c     An integer specifying the number of CUPs used. Default is 8
        -h     Show help message
'''

TEMP_Insertion = {
    'x' : 30,
    'm' : 1
}