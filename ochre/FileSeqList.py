from oseq.Sequence import Seq
from oseq.SeqList import SeqList
from oseq.FileFormats import FASTA, FASTQ, guess_filetype
import gzip
import bz2
import os
import itertools


class FileSeqList(SeqList):
    def __init__(self, fh, file_type=None, qual_type='guess'):
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

    def __len__(self):
        #TODO: is there a faster way?
        return sum(1 for _ in self)

    def __iter__(self):
        enc = 'utf-8'
        self._file.seek(0)
        if self._ftype == 'fa':
            file_reader = FASTA.read(self._file, enc)
        elif self._ftype == 'fq':
            #TODO: should be able to direct the type of quality score
            file_reader = FASTQ.read(self._file, enc, self._qtype)
        for seq in file_reader:
            yield seq

    def get_filename(self, frmt=None, tmp_dir='/tmp/ochre/'):
        if frmt is None:
            return os.path.abspath(self._file.name)
        else:
            if isinstance(frmt, str):
                frmt == (frmt,)
            if self._ftype in frmt:
                return os.path.abspath(self._file.name)
            else:
                file_name = os.path.join(tmp_dir, \
                  'seqs-' + str(id(self)) + '.' + frmt[0])
                self.write(file_name, frmt[0])
                return file_name

    def properties(self):
        pass


class PairedFileSeqList(object):
    def __init__(self, fh1, fh2, file_type=None, qual_type='guess', \
      seq_length=0):
        self._f1 = FileSeqList(fh1, file_type=file_type, qual_type=qual_type)
        self._f2 = FileSeqList(fh2, file_type=file_type, qual_type=qual_type)
        if self._f1._ftype != self._f2._ftype:
            raise TypeError('Files 1 and 2 are not the same type.')

        self.seq_length = seq_length

    def __len__(self):
        return sum(1 for _ in self)

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
        """Returns the two sets of sequences interleaved with each other
           for input to Velvet and other programs."""
        for s1, s2 in zip(self._f1, self._f2):
            yield s1
            yield s2
