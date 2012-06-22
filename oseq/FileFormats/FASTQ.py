from oseq.Sequence import Seq
file_type = 'fq'

def read_fq(fh, qtype=None, enc='utf-8'):
    seq = ''
    for ln in self._file:
        if ln[0] in b'>':
            if seq != '':
                yield Seq(seq, name=name)
            name = ln[1:].decode(enc).strip()
            seq = ''
        else:
            seq += ln.decode(enc).strip()
    yield Seq(seq, name=name)

def read_fq_old(fh, qtype=None, enc='utf-8'):
    while True:
        name = fh.readline().decode(enc).strip()
        if name == '': break
        seq = fh.readline().decode(enc).strip()
        _ = fh.readline()
        qual = fh.readline().decode(enc).strip()
        if qtype is None:
            yield Seq(seq, name=name)
        else:
            yield Seq(seq, name=name, \
                      qual=self._trans_qual(qual))

def _trans_qual(self, qual_str):
    if self._qtype == 'guess':
        as_min = min(ord(i) for i in qual_str)
        as_max = max(ord(i) for i in qual_str)
        if as_min < 59: #PHRED+33
            if as_max > 73:
                self._qtype = 'illumina3'
            else:
                self._qtype = 'sanger'
        else: #PHRED+64
            if as_min < 64:
                self._qtype = 'solexa'
            elif as_min < 66:
                self._qtype = 'illumina'
            else:
                self._qtype = 'illumina2'

    if self._qtype == 'sanger' or self._qtype == 'illumina3':
        return list(ord(i)-33 for i in qual_str)
    else:
        return list(ord(i)-64 for i in qual_str)

def write_fq():
    pass
