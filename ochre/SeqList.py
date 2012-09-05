import os
import io
import itertools
import types
from ochre.Sequence import Seq


class SeqList(object):
    """
    A Python object that represents a collection of sequences.
    The underlying storage can be in a file (through the
    subclass FileSeqList) or in a list of Seq objects.
    """
    def __new__(cls, sqlst, *args, **kwargs):
        from ochre.FileSeqList import FileSeqList
        if isinstance(sqlst, str):
            if os.path.extsep in sqlst or os.path.sep in sqlst:  # probably a file
                #TODO: check if URL
                _file = open(sqlst, 'r')
                return super(SeqList, cls).__new__(FileSeqList, _file, *args, **kwargs)
        elif isinstance(sqlst, io.IOBase):
            return super(SeqList, cls).__new__(FileSeqList, sqlst, *args, **kwargs)
        # for python 2, we need to do the following check
        try:
            if isinstance(sqlst, file):
                return super(SeqList, cls).__new__(FileSeqList, sqlst, *args, **kwargs)
        except NameError:
            pass
        # end python 2 code
        return super(SeqList, cls).__new__(SeqList, sqlst, *args, **kwargs)

    def __init__(self, sqlst, loose_indexing=False, other=None):
        """
        Create a list of sequences.

        Parameters:
            sqlst (file, str, list): data to create the list of
                sequences from. If *sqlst* is a file or string
                containing a path delimiter (typically / or .)
                read that file and return the sequences within.
                If *sqlst* is a string, separate it by comma
                and return each sequence. If *sqlst* is a list,
                make each item into a sequence.
            loose_indexing (bool): if true, match sequences that
                start with a given string and not the entire
                string, i.e. seqlist['te'] will match a sequence
                with the key, 'test'.
            other_csv (file, str): a csv file containing a list
                of data associated with each sequence, e.g.
                coverages, OTUs, etc. in the format:
                    "SeqName, Data1, Data2, ..."
                Data1, Data2, etc can then be obtained from the
                sequence through the info parameter.
        """
        if isinstance(sqlst, str):
            self._seqs = [Seq(s) for s in sqlst.split(',')]
        elif isinstance(sqlst, list):
            if isinstance(sqlst[0], str):
                self._seqs = [Seq(s) for s in sqlst]
            elif isinstance(sqlst[0], Seq):
                self._seqs = sqlst
        elif isinstance(sqlst, itertools.islice) or \
          isinstance(sqlst, types.GeneratorType):
            self._seqs = list(sqlst)
        elif sqlst is not None:
            # sqlst should be none if this is getting called
            # from FileSeqList, because that's already handling
            # the same things with _file
            raise TypeError('sqlst is not a file, str, or list.')

        if other is not None:
            if isinstance(sqlst, str):
                other = open(other, 'r')
            keys = other.readline().split(',')[1:]
            self._other = {}
            for ln in other:
                sep = ln.split(',')[1:]
                self._other[1:] = dict(zip(keys, sep[1:]))
            if self._seqs is not None:
                for s in self._seqs:
                    s.info = self._other.get(s.name, {})
        else:
            self._other = None

        self.loose_indexing = loose_indexing

    def __len__(self):
        return len(self._seqs)

    def __iter__(self):
        for s in self._seqs:
            yield s

    def __getitem__(self, key):
        if isinstance(key, str):
            for itm in self:
                if itm.name == key:
                    return itm
                elif itm.name.startswith(key) and self.loose_indexing:
                    return itm
            raise KeyError(str(key))
        elif isinstance(key, slice):
            return SeqList(itertools.islice(iter(self), \
                           key.start, key.stop, key.step))
        else:
            try:
                return next(itertools.islice(iter(self), key, key + 1))
            except StopIteration:
                raise IndexError("List index out of range.")

    def get_file(self, frmt='fa', fdir=None):
        """
        Write my sequences to disk in a file *seqs-######* in the
        directory *fdir* making sure they're in a format specified by the list "frmt"
        and return a FileSeqList.

        Useful for quickly getting a file in the correct format to
        pass onto an external program.

        Parameters:
            frmt (str): a file type (e.g. FASTA, FASTQ, etc)
                (default: FASTA)
            fdir (str): a directory to save the file in
                (default: the temp directory)
        Returns:
            the file name of the saved sequences
        """
        from ochre.External.External import temp

        if isinstance(frmt, str):
            frmt = (frmt,)

        fname = 'seqs-' + str(id(self)) + '.' + frmt[0]
        if fdir is None:
            file_name = temp(fname)
        else:
            file_name = os.path.join(fdir, fname)

        self.write(file_name, frmt[0])
        return file_name

    def write(self, filename, file_type=None):
        """
        Write out my sequences to a file.

        Parameters:
            filename (str): the name of the file to write to
            file_type (str): a file type (e.g. FASTA, FASTQ, etc)
                (default: determine from the file name,
                 otherwise FASTA)
        """
        from ochre.FileSeqList import FileSeqList
        from ochre.FileFormats.FileFormats import guess_filetype, file_writer

        if hasattr(filename, 'write'):
            fh = filename
            filename = fh.name
        else:
            fh = open(filename, 'w')

        if file_type is None:
            ftype = guess_filetype(filename.split(os.path.extsep)[-1])
            if ftype == '':
                ftype = 'fa'  # return a FASTA file if nothing else
        else:
            ftype = guess_filetype(file_type)

        file_writer(ftype, fh, self)

        if filename != '<stdout>':
            fh.close()
            return FileSeqList(filename)

    def n(self, interval=0.5):
        """
        Calculate N-values (i.e. N50) from a group of sequences

        Parameters:
            interval (float): a percentage
        Returns:
            *interval* percentage of the bases in this file are
            in a sequence of this length or longer
        """
        seq_len = [len(seq) for seq in self]
        seq_len.sort(reverse=True)
        s = sum(seq_len)
        limit = s * (1 - interval)
        for l in seq_len:
            s -= l
            if s <= limit:
                return l

    def filter(self, min_len=None, max_len=None, min_gc=None, \
      max_gc=None, min_cov=None, max_cov=None):
        pass

    def sample(self, n, replace=False):
        """
        Randomly sample sequences from this seqlist.

        Parameters:
            n (int): number of sequences to sample
            replace (bool): sample with replacement?
                (default: False)
        Returns:
            SeqList with *n* sampled sequences.
        """
        import random

        if replace:
            return SeqList(random.choice(list(self)) for _ in range(n))
        else:
            return SeqList(random.sample(list(self), n))

    def stats(self):
        """'Quick' statistics on sequences."""
        seqs, bps = len(self), sum(len(s) for s in self)
        o = 'Sequences: ' + str(seqs) + '\n'
        o += 'Basepairs: ' + str(bps) + '\n\n'
        o += 'Shortest: ' + str(min(len(s) for s in self)) + ' bp\n'
        o += 'N90: ' + str(self.n(0.9)) + ' bp\n'
        o += 'Average: ' + str(bps / seqs) + ' bp\n'
        o += 'N50: ' + str(self.n(0.5)) + ' bp\n'
        o += 'N10: ' + str(self.n(0.1)) + ' bp\n'
        o += 'Longest: ' + str(max(len(s) for s in self)) + ' bp'
        return o
