from ochre.FileFormats.FileFormats import SeqFile


class FASTQ(SeqFile):
    file_type = 'fq'
    abbrevs = ('fastq', 'fq')

    #def read_fq(fh, qtype=None, enc='utf-8'):
    #    seq, qual = '',''
    #    for ln in self._file:
    #        if ln[0] == b'@':
    #            if seq != '':
    #                yield Seq(seq, name=name)
    #            name = ln[1:].decode(enc).strip()
    #            seq = ''
    #        elif ln[0] == b'+':
    #        else:
    #            seq += ln.decode(enc).strip()
    #    yield Seq(seq, name=name)

    def read(self, fh, enc='utf-8', qtype=None):
        fh.seek(0)
        while True:
            name = fh.readline().decode(enc).strip()
            if name == '':
                break
            seq = fh.readline().decode(enc).strip()
            fh.readline()
            if qtype is None:
                fh.readline()  # quality string, can skip
                yield seq, name, {}
            else:
                qual_str = fh.readline().decode(enc).strip()
                if qtype == 'guess':
                    qtype = self._guess_qual(qual_str)
                if qtype in ['sanger', 'illumina3']:
                    qual = list(ord(i) - 33 for i in qual_str)
                else:
                    qual = list(ord(i) - 64 for i in qual_str)
                yield seq, name, {'qual': qual}


    def _guess_qual(self, qual_str):
        as_min = min(ord(i) for i in qual_str)
        as_max = max(ord(i) for i in qual_str)
        if as_min < 59:  # PHRED+33
            if as_max > 73:
                qtype = 'illumina3'
            else:
                qtype = 'sanger'
        else:  # PHRED+64
            if as_min < 64:
                qtype = 'solexa'
            elif as_min < 66:
                qtype = 'illumina'
            else:
                qtype = 'illumina2'
        return qtype


    def write(self, fh, seqs, qtype=None):
        for seq in seqs:
            fh.write('@' + seq.name + '\n')
            fh.write(seq.seq + '\n')
            fh.write('+\n')
            if seq.qual is None:
                #make up a quality score
                qual_str = len(seq.seq) * chr(35)
            elif qtype in ['sanger', 'illumina3']:
                qual_str = ''.join(chr(i + 33) for i in seq.qual)
            else:
                qual_str = ''.join(chr(i + 64) for i in seq.qual)
            fh.write(qual_str + '\n')
