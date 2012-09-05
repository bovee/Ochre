#Trimmomatic(, leading=3,trailing=3,slidingwin=(4,15),minlen=36)

def Trimmomatic(seqs, leading=None, trailing=None, slidingwin=None, minlen=None):
    c = [app('JAVA','java'),'-classpath',app('TRIMMOMATIC','trimmomatic-0.20.jar'),'org.usadellab.trimmomatic.TrimmomaticPE','-phred33']

    """${DATA_DIR}/MHL_7M_R1.fastq.gz
    ${DATA_DIR}/MHL_7M_R2.fastq.gz
    ${DATA_DIR}/MHL_7M_PAIRED_1.fastq.gz
    ${DATA_DIR}/MHL_7M_UNPAIRED_1.fastq.gz
    ${DATA_DIR}/MHL_7M_PAIRED_2.fastq.gz
    ${DATA_DIR}/MHL_7M_UNPAIRED_2.fastq.gz"""

    """ILLUMINACLIP:${DATA_DIR}/illumina.fa:2:40:15"""
    if leading is not None:
        c += ['LEADING:'+str(leading)]
    if trailing is not None:
        c += ['TRAILING:'+str(trailing)]
    if slidingwin is not None:
        c += ['SLIDINGWINDOW:'+str(slidingwin[0])+':'+str(slidingwin[1])]
    if headcrop is not None:
        c += ['HEADCROP:'+str(headcrop)]
    if minlen is not None:
        c += ['MINLEN:'+str(minlen)]
    run([c])
