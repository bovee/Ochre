file_type = 'fa'


def read(fh, enc='utf-8'):
    name, seq = '', ''
    for ln in fh:
        if ln[0] in b'>':
            if seq != '':
                yield seq, name, {}
            name = ln[1:].decode(enc).strip()
            seq = ''
        else:
            seq += ln.decode(enc).strip()
    yield seq, name, {}


def write(fh, seqs):
    for seq in seqs:
        fh.write('>' + seq.name + '\n')
        fh.write(seq.seq + '\n')
