from oseq.Sequence import Seq
from oseq.SeqList import SeqList
from oseq.FileFormats import FASTA, FASTQ, guess_filetype
import gzip
import bz2
import os

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
        filetype = guess_filetype(ext,mgc[0])
        if filetype == '':
            pass
        elif filetype in ['em','gb']:
            raise NotImplementedError
        self._ftype = filetype

    #def __len__(self):
        #TODO: figure out what to do here
        #return sum(1 for i in iter(self))
    #    return 0

    def __iter__(self):
        enc = 'utf-8'
        self._file.seek(0)
        if self._ftype == 'fa':
            file_reader = FASTA.read_fa(self._file, enc)
        elif self._ftype == 'fq':
            #TODO: should be able to direct the type of quality score
            file_reader = FASTQ.read_fq(self._file, enc, self._qtype)
        for seq in file_reader:
            yield seq

    def get_filename(self, frmt=None, tmp_dir='/tmp/ochre/'):
        if frmt is None:
            return os.path.abspath(self._file.name)
        else:
            if isinstance(frmt,str):
                frmt == (frmt,)
            if self._ftype in frmt:
                return os.path.abspath(self._file.name)
            else:
                file_name = os.path.join(tmp_dir, \
                                         'seqs-'+str(id(self))+'.'+frmt[0])
                self.write(file_name, frmt[0])
                return file_name

    def properties(self):
        pass
