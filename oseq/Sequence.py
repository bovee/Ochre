from oseq.Misc import tTables

class Seq(object):
    """A representation of the sequence of a biological macromolecule."""
    def __new__(cls, seq, *args, seq_type=None, **kwargs):
        if cls is Seq: #not a subclass, so guess which type it is
            if (set(seq.upper()).issubset(set('AUTGCN-')) and \
              seq_type is None) or seq_type == 'DNA' or seq_type == 'RNA':
                from oseq.NucSeq import NASeq
                return super(Seq, cls).__new__(NASeq, seq, *args, \
                                               seq_type=seq_type, **kwargs)
            elif seq_type == 'PROTEIN':
                from oseq.PepSeq import PepSeq
                return super(Seq, cls).__new__(PepSeq, seq, *args, \
                                               seq_type=seq_type, **kwargs)

        return super(Seq,cls).__new__(Seq, seq, *args, **kwargs)

    def __init__(self, seq, *args, name='', seq_type=None, \
                     qual=None, case_qual=False, **kwargs):
        self.seq = seq.upper()
        self.qual = qual
        self.name = name

        #the case of the letters indicates the quality of the sequence
        if case_qual:
            self.qual = [40 if i.isupper() else 20 for i in seq]

        #sets the type of the sequence if given, otherwise figures it out later
        self._stype = seq_type

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
    def reverse(self, **kwargs):
        """Returns the sequence reversed."""
        return Seq(self.seq[::-1])

    def translate(self, to_type='PROTEIN', from_type=None, table='standard'):
        """Change a sequence into another sequence."""
        assert to_type in ['DNA','RNA','PROTEIN']

        if from_type is None:
            from_type = self.stype
        if from_type == to_type:
            return self

        if to_type == 'PROTEIN':
            tab = tTables(table)
            flatten = lambda l: [item for sublist in l for item in sublist]
            rtab = dict(flatten([zip(tab[i],len(tab[i])*[i]) for i in tab]))

            if 'RNA' in from_type:
                rtab = dict([(i.replace('T','U'),rtab[i]) for i in rtab])

            def codons():
                for i in range(0,len(self.seq)-len(self.seq)%3,3):
                    yield self.seq[i:i+3]

            rseq = ''.join([rtab.get(cdn, 'X') for cdn in codons()])

            return Seq(rseq, seq_type='PROTEIN')
        else: #translate to NA
            if from_type == 'PROTEIN':
                tab = tTables(table)
                try:
                    rseq = ''.join([tab[i][0] for i in self.seq])
                except KeyError:
                    raise KeyError('Amino acid not found in table.')
                if to_type == 'RNA':
                    return Seq(rseq.replace('T','U'),seq_type='RNA')
                else:
                    return Seq(rseq,seq_type='DNA')
            elif from_type == 'DNA':
                #already handled the DNA->DNA case above
                return Seq(self.seq.replace('T','U'),seq_type='RNA')
            else:
                return Seq(self.seq.replace('U','T'),seq_type='DNA')

    def mass(self):
        raise NotImplementedError

    def freqs(self, lngth=1, rev_comp=True, overlapping=True):
        """ Calculate the relative abundance of all possible kmers.

        """
        import itertools

        sseqs = (self.seq[c:c+lngth] for c in range(len(self.seq)-lngth+1))
        ss_abun = {}

        if rev_comp:
            #all possible letter combinations
            pss = [''.join(i) for i in itertools.product(*(lngth*['ATGC']))]
            seq_map = {}
            for s in pss:
                if not self._invert(s[::-1]) in seq_map:
                    seq_map[s] = s
                    seq_map[invert(s[::-1])] = s

            #subsequences generator
            for ss in sseqs:
                if ss in ss_abun:
                    ss_abun[seq_map[ss]] += 1
                else:
                    ss_abun[seq_map[ss]] = 1
        else:
            #subsequences generator
            for ss in sseqs:
                if ss in ss_abun:
                    ss_abun[ss] += 1
                else:
                    ss_abun[ss] = 1
        return ss_abun
