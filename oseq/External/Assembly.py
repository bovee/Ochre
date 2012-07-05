def Velvet(seqs, kmer=25):
    if not (isinstance(seqs, tuple) or isinstance(seqs, list)):
        seqs = [seqs]
    for seq_file in seqs:
        pass
    c = [['velveth']]
    c += [['velvetg']]

def get_velvet_bam_file(velvet_afg):
    amos_bank = 'data.bnk'
    bout = 'data'
    c = [[app('AMOS','bank-transact'),'-m',velvet_afg,'-b',amos_bank,'-c']]
    c += [[app('AMOS','bank2fasta'),'-i','-b',amos_bank,'>',bout+'.fa']]
    c += [[app('AMOS','bank2contig'),'-i','-b',amos_bank,'>',bout+'.sam']]
    c += [[app('SAMTOOLS','samtools'),'faidx',bout+'.fa']]
    c += [[app('SAMTOOLS','samtools'),'import',bout+'.fa.fai', \
      bout+'.sam',bout+'.unsorted.bam']]
    c += [[app('SAMTOOLS','samtools'),'sort',bout+'.unsorted.bam',bout+'.bam']]
    c += [[app('SAMTOOLS','samtools'),'index',bout+'.bam']]
    run(c)
    return bout+'.bam'

def Newbler(seqs):
    pass
