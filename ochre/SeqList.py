from oseq.Sequence import Seq
import os
import io
import itertools


class SeqList(object):
    def __new__(cls, sqlst, *args, **kwargs):
        from oseq.FileSeqList import FileSeqList
        if isinstance(sqlst, str):
            if os.path.extsep in sqlst or os.path.sep in sqlst:  # probably a file
                #TODO: check if URL
                _file = open(sqlst, 'r')
                return super(SeqList, cls).__new__(FileSeqList, _file, *args, **kwargs)
        elif isinstance(sqlst, io.IOBase):
            return super(SeqList, cls).__new__(FileSeqList, sqlst, *args, **kwargs)
        return super(SeqList, cls).__new__(SeqList, sqlst, *args, **kwargs)

    def __init__(self, sqlst):
        if isinstance(sqlst, str):
            self._seqs = [Seq(s) for s in sqlst.split(',')]
        elif isinstance(sqlst, list):
            self._seqs = [Seq(s) for s in sqlst]
        elif isinstance(sqlst, itertools.islice):
            self._seqs = list(sqlst)

    def __len__(self):
        return len(self._seqs)

    def __iter__(self):
        for s in self._seqs:
            yield s

    def __getitem__(self, key):
        if isinstance(key, str):
            raise NotImplementedError
        elif isinstance(key, slice):
            return SeqList(itertools.islice(iter(self), \
                           key.start, key.stop, key.step))
        else:
            try:
                return next(itertools.islice(iter(self), key, key + 1))
            except StopIteration:
                raise IndexError("list index out of range")

    def get_file(self, frmt='fa', dir=None):
        from oseq.FileSeqList import FileSeqList
        if isinstance(frmt, str):
            frmt == (frmt,)
        if isinstance(self, FileSeqList):
            if self._ftype in frmt and dir is None:
                return self._file.name
            else:
                dir = '/tmp/ochre'
        else:
            dir = '/tmp/ochre'
        file_name = os.path.join(dir, 'seqs-' + str(id(self)) + '.' + frmt[0])
        self.write(file_name, frmt[0])
        return file_name

    def write(self, filename, file_type=None):
        if file_type is None:
            file_type = filename.split(os.path.extsep)[-1]

        fh = open(filename, 'w')
        if file_type == 'fasta' or file_type == 'fa':
            from oseq.FileFormats import FASTA
            FASTA.write(fh, self._seqs)
        elif file_type == 'fastq' or file_type == 'fq':
            from oseq.FileFormats import FASTQ
            FASTQ.write(fh, self._seqs)
        self._file = filename
        self._ftype = file_type
        #TODO: make me into a FileSeqList now?

    def n(self, interval):
        """Calculate N-values (i.e. N50) from a group of sequences"""
        seq_len = [len(seq) for seq in self]
        seq_len.sort(reverse=True)
        s = sum(seq_len)
        limit = s * (1 - interval)
        for l in seq_len:
            s -= l
            if s <= limit:
                return l

    def stats(self):
        """'Quick' statistics on sequences in here."""
        seqs, bps = len(self), sum(len(s) for s in self)
        print('Sequences: ' + str(seqs))
        print('Basepairs: ' + str(bps))
        print()
        print('Shortest: ' + str(min(len(s) for s in self)) + ' bp')
        print('N90: ' + str(self.n(0.9)) + ' bp')
        print('Average: ' + str(bps / seqs) + ' bp')
        print('N50: ' + str(self.n(0.5)) + ' bp')
        print('N10: ' + str(self.n(0.1)) + ' bp')
        print('Longest: ' + str(max(len(s) for s in self)) + ' bp')
