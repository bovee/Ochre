from ochre.SeqList import SeqList
from ochre.External.External import app, temp, run
import os.path as op


def BLAST(seq, seqs, evalue=10):
    lib_dir = temp('blast-' + str(id(seqs)))
    res_file = temp('blast-' + str(id(seqs)), str(id(seq)))
    dbtype = 'nucl'  # prot

    if op.exists(op.join(lib_dir, seqs.get_file('fa') + '.phr')) \
      or op.exists(op.join(lib_dir, seqs.get_file('fa') + '.nhr')):
        pass

    c = [[app('BLAST','makeblastdb'),'-in',seqs.get_file('fa'),'-out',lib_dir,'-parse_seqids','-dbtype',dbtype]]
    c += [[app('BLAST','blastn'),'-db',lib_dir,'-query',seq.get_file('fa'),'-outfmt','"6 sseqid"','-out',res_file+'-out','-evalue',str(evalue)]]
    c += [[app('BLAST','blastdbcmd'),'-db',lib_dir,'-entry_batch',res_file+'-out','-out',res_file+'.fa']]
    run(c)
    return SeqList(open(res_file + '.fa', 'r'))
