import io
import Bio.SeqIO
import Bio.AlignIO
from ochre.FileFormats.FileFormats import SeqFile


class GenBank(SeqFile):
    file_type = 'gb'
    abbrevs = ('gb', 'genbank')

    def read(self, fh, enc='utf-8', *args, **kwargs):
        fh.seek(0)
        fth = io.TextIOWrapper(fh)
        for rcd in Bio.SeqIO.parse(fth, 'genbank'):
            yield str(rcd.seq), rcd.id, {}

    def write(self, fh, seqs):
        raise NotImplementedError


class Clustal(SeqFile):
    file_type = 'ctl'
    abbrevs = ('aln', 'ctl', 'clustal')

    def read(self, fh, enc='utf-8', *args, **kwargs):
        fh.seek(0)
        fth = io.TextIOWrapper(fh)
        for msa in Bio.AlignIO.parse(fth, 'clustal'):
            for rcd in msa:
                yield str(rcd.seq), rcd.id, {}

    def write(self, fh, seqs):
        raise NotImplementedError


class Stockholm(SeqFile):
    file_type = 'sto'
    abbrevs = ('sto', 'stk', 'stockholm')

    def read(self, fh, enc='utf-8', *args, **kwargs):
        fh.seek(0)
        fth = io.TextIOWrapper(fh)
        for rcd in Bio.AlignIO.parse(fth, 'stockholm'):
            yield str(rcd.seq), rcd.id, {}

    def write(self, fh, seqs):
        raise NotImplementedError
