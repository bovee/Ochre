from ochre.External.External import app, temp, run
import subprocess
import re
import os.path as op
from ochre.FileSeqList import FileSeqList, PairedFileSeqList


def Velvet(mseqs, kmer=25, out_dir=None):
    if not (isinstance(mseqs, tuple) or isinstance(mseqs, list)):
        mseqs = [mseqs]

    a = subprocess.getoutput(app('VELVET', 'velveth'))
    max_kmer_len = re.search('MAXKMERLENGTH.+?(\d+)\n', a).group(1)
    if kmer > max_kmer_len:
        raise ValueError('K-mer value too high for Velvet. Recompile Velvet.')
    #categories = re.search('CATEGORIES.+?(\d+)\n', a).group(1)

    if out_dir is None:
        out_dir = temp('velvet-k' + str(kmer) + '-' + \
          '-'.join(str(id(seq_file)) for seq_file in mseqs))

    for seq_file in mseqs:
        if sum(len(s) for s in seq_file[:10]) / 10.0 > 150:
            opt = '-long'  # these are 454 or longer reads
        else:
            opt = '-short'

        if isinstance(seq_file, PairedFileSeqList):
            opt += 'Paired'

    c = [[app('VELVET', 'velveth'), out_dir, '', '']]
    c += [[app('VELVET', 'velvetg'), out_dir]]
    run(c)
    return FileSeqList(op.join(out_dir, 'contigs.fa'))


def get_velvet_bam_file(velvet_afg):
    #TODO: this might not work?
    amos_bank = 'data.bnk'
    bout = 'data'
    c = [[app('AMOS', 'bank-transact'), '-m', velvet_afg, '-b', amos_bank, '-c']]
    c += [[app('AMOS', 'bank2fasta'), '-i', '-b', amos_bank, '>', bout + '.fa']]
    c += [[app('AMOS', 'bank2contig'), '-i', '-b', amos_bank, '>', bout + '.sam']]
    c += [[app('SAMTOOLS', 'samtools'), 'faidx', bout + '.fa']]
    c += [[app('SAMTOOLS', 'samtools'), 'import', bout + '.fa.fai', \
      bout + '.sam', bout + '.unsorted.bam']]
    c += [[app('SAMTOOLS', 'samtools'), 'sort', bout + '.unsorted.bam', bout + '.bam']]
    c += [[app('SAMTOOLS', 'samtools'), 'index', bout + '.bam']]
    run(c)
    return bout + '.bam'


def IDBA(seqs, out_dir=None):
    if out_dir is None:
        out_dir = temp('idba' + '-' + str(id(seqs)))
    fname = seqs.get_file('fa')
    c = [[app('IDBA_UD', 'idba_ud'), '-r', fname, '-o', out_dir]]
    # '--num_threads','8'
    run(c)


def Newbler(seqs, out_dir=None):
    if out_dir is None:
        out_dir = temp('newbler' + '-' + str(id(seqs)))
    fname = seqs.get_file('sff')
    c = [[app('NEWBLER', 'runAssembly'), '-o', out_dir, fname]]
    run(c)
    return FileSeqList(op.join(out_dir, '454AllContigs.fna'))
