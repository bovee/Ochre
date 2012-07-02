from oseq.FileSeqList import FileSeqList
from subprocess import call
import os.path as op

_tmp_dir = '/tmp/ochre'
_BLAST_dir = ''
_BWA_dir = ''
_TRIMMOMATIC_dir = ''
_VELVET_dir = ''

def Velvet(seqs, kmer=25, ins_length=300):
    if isinstance(seqs, tuple) or isinstance(seqs, list):
        pass
    call(['velveth'])
    call(['velvetg'])

def BWA(ref_seq, seqs):
    if True:
        call(['samse'])
    else:
        call(['sampe'])

def Trimmomatic(seqs, leading=3, trailing=3, slidingwin=(4,15), minlen=36):
    pass

def BLAST(seq, seqs, evalue=10):
    lib_dir = op.join(_tmp_dir,str(id(seqs)))
    res_file = op.join(lib_dir,str(id(seq)))
    dbtype = 'nucl' #prot

    c = [op.join(_BLAST_dir,'makeblastdb'),'-in',seqs.get_file('fa'),'-out',lib_dir,'-parse_seqids','-dbtype',dbtype]
    c += [op.join(_BLAST_dir,'blastn'),'-db',lib_dir,'-query',seqs.get_file('fa'),'-outfmt','"6 sseqid"','-out',res_file+'-out','-evalue',str(evalue)]
    c += [op.join(_BLAST_dir,'blastdbcmd'),'-db',lib_dir,'-entry_batch',res_file+'-out','-out',res_file+'.fa']
    run_cmds(c)
    return FileSeqList(open(res_file+'.fa','r'))
