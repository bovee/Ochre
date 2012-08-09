from ochre.Sequence import Seq
from ochre.SeqList import SeqList
from ochre.FileFormats import guess_filetype, file_reader
import gzip
import bz2
import os
import itertools


class FileSeqList(SeqList):
    def __init__(self, fh, file_type=None, qual_type='guess', loose_indexing=False):
        if isinstance(fh, str):
            fh = open(fh, 'rb')

        if 'b' not in fh.mode:
            raise IOError("File not opened as binary.")

        self._qtype = qual_type
        ext = fh.name.split(os.path.extsep)[-1]
        mgc = fh.read(2)
        fh.seek(0)

        #transparently decompress the file
        if ext == 'gz' or ext == 'gzip' or mgc == b'\x1f\x8b':
            self._file = gzip.GzipFile(fileobj=fh)
            ext = self._file.name.split(os.path.extsep)[-2]
            mgc = self._file.read(2)
            self._file.seek(0)
        elif ext == 'bz' or mgc == b'\x42\x5a':
            raise NotImplementedError
        else:
            self._file = fh

        #now deal with how to handle the actual data
        filetype = guess_filetype(ext, mgc[0])
        if filetype == '':
            pass
        elif filetype in ['em', 'gb']:
            raise NotImplementedError
        self._ftype = filetype

        self.loose_indexing = loose_indexing

    def _raw_reads(self):
        enc = 'utf-8'
        self._file.seek(0)
        return file_reader(self._ftype, self._file, enc, self._qtype)

    def __len__(self):
        return sum(1 for _ in self._raw_reads())

    def __iter__(self):
        for seq, name, other in self._raw_reads():
            if 'qual' in other.keys():
                yield Seq(seq, name=name, qual=other['qual'])
            else:
                yield Seq(seq, name=name)

    def get_file(self, frmt='fa', fdir=None):
        if isinstance(frmt, str):
            frmt == (frmt,)

        if self._ftype in frmt and fdir is None:
            return self  # os.path.abspath(self._file.name)
        #elif self._ftype in frmt:  # copy the file to the right place
        #    pass
        else:
            return super(FileSeqList, self).get_file(frmt, fdir)

    def properties(self):
        pass


class PairedFileSeqList(object):
    def __init__(self, fh1, fh2, file_type=None, qual_type='guess', \
      seq_length=0):
        """
        Parameters:
            fh1 (str, file): the forward reads
            fh2 (str, file): the backward reads
            file_type (str): file format (e.g. FASTA, FASTQ, etc)
            qual_type (str): quality type (e.g. PHRED+64, etc)
            seq_length (int): length of each sequence

        Returns:
            a PairedFileSeqList object
        """
        self._f1 = FileSeqList(fh1, file_type=file_type, qual_type=qual_type)
        self._f2 = FileSeqList(fh2, file_type=file_type, qual_type=qual_type)
        if self._f1._ftype != self._f2._ftype:
            raise TypeError('Files 1 and 2 are not the same type.')

        self.seq_length = seq_length

    def __len__(self):
        return len(self._f1) + len(self._f2)

    def __iter__(self):
        for s1, s2 in zip(self._f1, self._f2):
            #calculate length of gap in middle
            gl = max(0, self.seq_length - len(s1.seq) - len(s2.seq))
            yield Seq(s1.seq + (gl * 'N') + s2.reverse(compliment=True).seq, \
              name=s1.name, qual=s1.qual + (gl * [0]) + s2.qual)

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

    def interleaved(self):
        """
        Return the two sets of sequences interleaved with each other
        for input to Velvet and other programs.
        """
        for s1, s2 in zip(self._f1, self._f2):
            yield s1
            yield s2
