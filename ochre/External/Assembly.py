from ochre.External.External import app, temp, run
import subprocess
import re
from ochre.FileSeqList import PairedFileSeqList


def Velvet(seqs=None, kmer=25):
    if not (isinstance(seqs, tuple) or isinstance(seqs, list)):
        seqs = [seqs]

    a = subprocess.getoutput(app('VELVET','velveth'))
    max_kmer_len = re.search('MAXKMERLENGTH.+?(\d+)\n', a).group(1)
    categories = re.search('CATEGORIES.+?(\d+)\n', a).group(1)

    tmp_dir = temp('velvet-k' + str(kmer) + '-' + \
      '-'.join(str(id(seq_file)) for seq_file in seqs))

    for seq_file in seqs:
        if sum(len(s) for s in seq_file[:10]) / 10.0 > 300:
            opt = 'long'  # these are 454 or longer reads
        else:
            opt = 'short'

        if isinstance(seq_file, PairedFileSeqList):
            opt += 'Paired'

    c = [[app('VELVET','velveth'),tmp_dir,'','']]
    c += [[app('VELVET','velvetg'),tmp_dir]]
    run(c)


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
