from oseq.FileSeqList import FileSeqList
import subprocess

_tmp_dir = '/tmp/'

def Velvet(seqs, paired_seqs=None, kmer=25, ins_length=300):
    pass

def BWA(ref_seq, seqs):
    pass

def Trimmomatic(seqs, paired_seqs=None, leading=3, trailing=3, slidingwin=(4,15), minlen=36):
    pass

def BLAST(seq, seqs, evalue=10):
    blast_dir = ''
    blast_results = ''
    dbtype = 'nucl'
    res_file = ''
    subprocess.call(['makeblastdb','-in',get_seqs_file(seqs,'fa'),'-out',blast_dir,'-parse_seqids','-dbtype',dbtype])
    subprocess.call(['blastn','-db',blast_dir,'-query',get_seqs_file(seqs,'fa'),'-outfmt','"6 sseqid"','-out',blast_results,'-evalue',str(evalue)])
    subprocess.call(['blastdbcmd','-db',blast_dir,'-entry_batch',blast_results,'-out',res_file])

def get_seqs_file(seqs, tmp_dir='/tmp/ochre', frmt='fa'):
    assert isinstance(seqs, SeqList)
    if isinstance(frmt, str):
        frmt == (frmt,)
    if isinstance(seqs, FileSeqList):
        if seqs._ftype in frmt:
            return file_name
    seqs.write(os.path.join(tmp_dir,'seqs.'+frmt[0]), frmt[0])
    return file_name
