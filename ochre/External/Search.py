from ochre.SeqList import SeqList
from ochre.NucSeq import NASeq
from ochre.External.External import app, temp, run
import os.path as op


def BLAST(seq, seqs, evalue=10):
    lib_dir = temp('blast-' + str(id(seqs)))
    res_file = temp('blast-' + str(id(seqs)), str(id(seq)))

    dbtype = 'nucl' if isinstance(seqs[0], NASeq) else 'prot'
    seqtype = 'nucl' if isinstance(seq, NASeq) else 'prot'

    if not op.exists(op.join(lib_dir, seqs.get_file('fa') + '.phr')) \
      and not op.exists(op.join(lib_dir, seqs.get_file('fa') + '.nhr')):
        c = [[app('BLAST','makeblastdb'),'-in',seqs.get_file('fa'),'-out',lib_dir,'-parse_seqids','-dbtype',dbtype]]
    else:
        c = []

    if dbtype == 'nucl' and seqtype == 'nucl':
        blprog = 'blastn'  # TBLASTX
    elif dbtype == 'nucl' and seqtype == 'prot':
        blprog = 'blastx'
    elif dbtype == 'prot' and seqtype == 'nucl':
        blprog = 'tblastn'
    elif dbtype == 'prot' and seqtype == 'prot':
        blprog = 'blastp'  # PSIBLAST, PHIBLAST

    c += [[app('BLAST',blprog),'-db',lib_dir,'-query',seq.get_file('fa'),'-outfmt','"6 sseqid"','-out',res_file+'-out','-evalue',str(evalue)]]
    c += [[app('BLAST','blastdbcmd'),'-db',lib_dir,'-entry_batch',res_file+'-out','-out',res_file+'.fa']]
    run(c)

    return SeqList(open(res_file + '.fa', 'r'))
