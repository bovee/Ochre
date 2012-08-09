import os
import io
import itertools
from ochre.Sequence import Seq


class SeqList(object):
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

    def __init__(self, sqlst, loose_indexing=False):
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
        """
        if isinstance(sqlst, str):
            self._seqs = [Seq(s) for s in sqlst.split(',')]
        elif isinstance(sqlst, list):
            if isinstance(sqlst[0], str):
                self._seqs = [Seq(s) for s in sqlst]
            elif isinstance(sqlst[0], Seq):
                self._seqs = sqlst
        elif isinstance(sqlst, itertools.islice):
            self._seqs = list(sqlst)

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
            a new FileSeqList representing the saved sequences
        """
        from ochre.FileSeqList import FileSeqList
        from ochre.External.External import temp

        if isinstance(frmt, str):
            frmt == (frmt,)

        fname = 'seqs-' + str(id(self)) + '.' + frmt[0]
        if fdir is None:
            file_name = temp(fname)
        else:
            file_name = os.path.join(fdir, fname)

        self.write(file_name, frmt[0])
        return FileSeqList(file_name)

    def write(self, filename, file_type=None):
        """
        Write out my sequences to a file.

        Parameters:
            filename (str): the name of the file to write to
            file_type (str): a file type (e.g. FASTA, FASTQ, etc)
                (default: determine from the file name)
        """
        if file_type is None:
            file_type = filename.split(os.path.extsep)[-1]

        fh = open(filename, 'w')
        if file_type == 'fasta' or file_type == 'fa':
            from ochre.FileFormats import FASTA
            FASTA.write(fh, self._seqs)
        elif file_type == 'fastq' or file_type == 'fq':
            from ochre.FileFormats import FASTQ
            FASTQ.write(fh, self._seqs)
        #TODO: make me into a FileSeqList now?
        #self._file = filename
        #self._ftype = file_type

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
