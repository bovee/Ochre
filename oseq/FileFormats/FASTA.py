from oseq.Sequence import Seq
file_type = 'fa'

def read_fa(fh, enc='utf-8'):
    seq = ''
    for ln in fh:
        if ln[0] in b'>':
            if seq != '':
                yield Seq(seq, name=name)
            name = ln[1:].decode(enc).strip()
            seq = ''
        else:
            seq += ln.decode(enc).strip()
    yield Seq(seq, name=name)

def write_fa(fh, seqs):
    for seq in seqs:
        fh.write('>'+seq.name+'\n')
        fh.write(seq.seq+'\n')

