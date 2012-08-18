class Seq(object):
    """A representation of the sequence of a biological macromolecule."""
    def __new__(cls, seq, *args, **kwargs):
        seq_type = kwargs.pop('seq_type', None)
        if cls is Seq:  # not a subclass, so guess which type it is
            if (set(seq.upper()).issubset(set('AUTGCN-')) and \
              seq_type is None) or seq_type == 'DNA' or seq_type == 'RNA':
                from ochre.NucSeq import NASeq
                return super(Seq, cls).__new__(NASeq, seq, *args, \
                                               seq_type=seq_type, **kwargs)
            elif seq_type == 'PROTEIN':
                from ochre.PepSeq import PepSeq
                return super(Seq, cls).__new__(PepSeq, seq, *args, \
                                               seq_type=seq_type, **kwargs)

        return super(Seq, cls).__new__(Seq, seq, *args, **kwargs)

    def __init__(self, seq, *args, **kwargs):
        """
        Parameters:
            seq (str):
            name (str):
            qual (list):
            seq_type (str):
            case_qual (bool):
        """

        self.seq = seq.upper()
        self.name = kwargs.get('name', '')
        self.qual = kwargs.get('qual', None)

        #the case of the letters indicates the quality of the sequence
        if kwargs.get('case_qual', False):
            self.qual = [40 if i.isupper() else 20 for i in seq]

        #sets the type of the sequence if given, otherwise figures it out later
        self._stype = kwargs.get('seq_type', None)

    @property
    def stype(self):
        if self._stype is not None:
            return self._stype
        elif set(self.seq).issubset(set('ATGC')):
            return 'DNA'
        elif set(self.seq).issubset(set('AUGC')):
            return 'RNA'
        elif set(self.seq).issubset(set('ATGC-')):
            return 'ALIGNED_DNA'
        elif set(self.seq).issubset(set('AUGC-')):
            return 'ALIGNED_RNA'
        else:
            return 'UNKNOWN'

    def __repr__(self):
        rstr = self.seq if len(self.seq) <= 56 else self.seq[:53] + '...'
        rstr += '' if self.qual is None else '!'
        if self.name != '':
            rstr += ' (' + (self.name if len(self.name) <= 20 \
                                else self.name[:17] + '...') + ')'
        return rstr

    def __str__(self):
        return self.seq

    def __len__(self):
        return len(self.seq)

    # Common Code
    def slid_win(self, size=4, overlapping=True):
        """Returns a sliding window along self.seq."""
        itr = iter(self.seq)
        if overlapping:
            buf = ''
            for _ in range(size):
                buf += next(itr)
            for l in itr:
                yield buf
                buf = buf[1:] + l
            yield buf
        else:
            buf = ''
            for l in itr:
                buf += l
                if len(buf) == size:
                    yield buf
                    buf = ''
            yield buf

    def reverse(self, **kwargs):
        """Returns the sequence reversed."""
        return Seq(self.seq[::-1])

    def mass(self):
        raise NotImplementedError

    def freqs(self, lngth=1, overlapping=True):
        """ Calculate the relative abundance of all possible kmers.

        """
        import itertools

        sseqs = (self.seq[c:c+lngth] for c in range(len(self.seq)-lngth+1))
        ss_abun = {}

        #subsequences generator
        for ss in sseqs:
            if ss in ss_abun:
                ss_abun[ss] += 1
            else:
                ss_abun[ss] = 1
        return ss_abun
