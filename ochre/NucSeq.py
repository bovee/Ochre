from ochre.Sequence import Seq


class NASeq(Seq):
    def gc(self):
        #TODO: fails if R or Y present?
        return float(sum([1 for i in self.seq if i in 'GC'])) / len(self.seq)

    def invert(self):
        return Seq(self.__invert__(self.seq, self._is_rna()))

    def reverse(self, compliment=False):
        """Returns the sequence reversed and possibly complimented."""
        if compliment:
            return Seq(self._invert(self.seq[::-1], self._is_rna()))
        else:
            return Seq(self.seq[::-1])

    def _is_rna(self):
        return True if 'RNA' in self.stype else False

    def _invert(self, seq, is_rna=False):
        invert_table = {'A': 'T', 'U': 'A', 'T': 'A', 'C': 'G', 'G': 'C',
                        'R': 'R', 'Y': 'Y', 'N': 'N', '-': '-'}
        if is_rna:
            invert_table['A'] = 'U'
        return ''.join(invert_table[i] for i in seq)

    def melting_temp(self, method='basic'):
        #http://bioinformatics.oxfordjournals.org/content/21/6/711.long
        import collections
        bps = collections.Counter(self.seq.upper())
        if method == 'basic':
            return 64.9 + 41.0 * (bps['G'] + bps['C'] - 16.4) / \
              (bps['A'] + bps['T'] + bps['G'] + bps['C'])

    def tetra_freqs(self, lgth=4, seq_map=None):
        import itertools
        #TODO: only good for DNA seqs

        def invert(seq):
            invert_table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            return ''.join([invert_table[i] for i in seq])

        def slid_win(itr, size=4):
            buf = ''
            for _ in range(size):
                buf += next(itr)
            for l in itr:
                yield buf
                buf = buf[1:] + l
            yield buf

        if seq_map is None:
            seq_map = {'': lgth * 'N'}
            for s in (''.join(i) for i in itertools.product(*(lgth * ['ATGC']))):
                if invert(s[::-1]) not in seq_map or s not in seq_map:
                    seq_map[s] = s
                    seq_map[invert(s[::-1])] = s

        #subsequences generator
        ss_abun = dict([(s, 0) for s in seq_map.values()])
        for ss in slid_win(iter(self.seq), lgth):
            ss_abun[seq_map.get(ss, lgth * 'N')] += 1
        return ss_abun
