from ochre.FileFormats.FileFormats import SeqFile


class FASTA(SeqFile):
    file_type = 'fa'
    abbrevs = ('fa', 'fasta', 'fna', 'faa')

    def read(self, fh, enc='utf-8', *args, **kwargs):
        fh.seek(0)
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

    def write(self, fh, seqs):
        for seq in seqs:
            fh.write('>' + seq.name + '\n')
            fh.write(seq.seq + '\n')
