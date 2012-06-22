from oseq.Sequence import Seq
import os
import io
import itertools
import inspect

class SeqList(object):
    def __new__(cls, sqlst, *args, **kwargs):
        from oseq.FileSeqList import FileSeqList
        if isinstance(sqlst, str):
            if os.path.extsep in sqlst or os.path.sep in sqlst: #probably a file
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
                return next(itertools.islice(iter(self),key,key+1))
            except StopIteration:
                raise IndexError("list index out of range")                
    
    def get_filename(self, frmt=None, tmp_dir='/tmp/ochre/'):
        if frmt is None:
            frmt = 'fa'
        elif isinstance(frmt, list) or isinstance(frmt, tuple):
            frmt = frmt[0]
        else:
            assert isinstance(frmt, str)
        file_name = os.path.join(tmp_dir,'seqs-'+str(id(self))+'.'+frmt[0])
        self.write(file_name, frmt[0])
        return file_name
    
    def properties(self):
        pass

    def align(self, method=''):
        '''Returns a SeqList of all of the Seqs in this list aligned.'''
        raise NotImplementedError

    def find(self, method=''):
        '''Searches through the SeqList and returns all similar Seqs.'''
        raise NotImplementedError

    def filter(self, method=''):
        raise NotImplementedError

    def assemble(self, method=''):
        raise NotImplementedError

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
