from oseq.SeqList import SeqList
from oseq.External.External import app, run_cmds
import os.path as op

tmp_dir = '/tmp/ochre'

def BLAST(seq, seqs, evalue=10):
    lib_dir = op.join(tmp_dir,str(id(seqs)))
    res_file = op.join(lib_dir,str(id(seq)))
    dbtype = 'nucl' #prot

    if op.exists(op.join(lib_dir,seqs.get_file('fa')+'.phr')):
        pass #only for prot sequences now? need *.nhr?


    c = [[app('BLAST','makeblastdb'),'-in',seqs.get_file('fa'),'-out',lib_dir,'-parse_seqids','-dbtype',dbtype]]
    c += [[app('BLAST','blastn'),'-db',lib_dir,'-query',seq.get_file('fa'),'-outfmt','"6 sseqid"','-out',res_file+'-out','-evalue',str(evalue)]]
    c += [[app('BLAST','blastdbcmd'),'-db',lib_dir,'-entry_batch',res_file+'-out','-out',res_file+'.fa']]
    run_cmds(c)
    return SeqList(open(res_file+'.fa','r'))
