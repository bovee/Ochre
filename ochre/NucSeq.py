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

    def to_nuc(self, table='standard', rna=False):
        if rna:
            return Seq(self.seq.replace('T', 'U'), seq_type='RNA')
        else:
            return Seq(self.seq.replace('U', 'T'), seq_type='DNA')

    def to_pep(self, table='standard'):
        from ochre.Misc import tTables
        tab = tTables(table)

        flatten = lambda l: [item for sublist in l for item in sublist]
        rtab = dict(flatten([zip(tab[i], len(tab[i]) * [i]) for i in tab]))

        if self._is_rna():
            rtab = dict([(i.replace('T', 'U'), rtab[i]) for i in rtab])

        rseq = ''.join(rtab.get(cdn, 'X') for cdn \
                in self.slid_win(3, overlapping=False) if len(cdn) == 3)

        return Seq(rseq, seq_type='PROTEIN')

    def melting_temp(self, method='basic'):
        #http://bioinformatics.oxfordjournals.org/content/21/6/711.long
        import collections
        bps = collections.Counter(self.seq.upper())
        if method == 'basic':
            return 64.9 + 41.0 * (bps['G'] + bps['C'] - 16.4) / \
              (bps['A'] + bps['T'] + bps['G'] + bps['C'])

    def nuc_freqs(self, lngth=4, seq_map=None):
        import itertools
        #TODO: only good for DNA seqs

        if seq_map is None:
            def invert(seq):
                invert_table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                return ''.join([invert_table[i] for i in seq])

            seq_map = {'': lngth * 'N'}
            for s in (''.join(i) for i in itertools.product(*(lngth * ['ATGC']))):
                if invert(s[::-1]) not in seq_map or s not in seq_map:
                    seq_map[s] = s
                    seq_map[invert(s[::-1])] = s

        #subsequences generator
        ss_abun = dict([(s, 0) for s in seq_map.values()])
        for ss in self.slid_win(lngth):
            ss_abun[seq_map.get(ss, lngth * 'N')] += 1
        return ss_abun

    def tetra_zscore(self, seq_map=None):
        import itertools
        from math import sqrt
        #TODO: only good for DNA seqs

        def rc(seq):  # reverse complement
            invert_table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
            return ''.join(invert_table.get(i, 'N') for i in seq[::-1])

        if seq_map is None:
            seq_map = {}
            seq_map[4] = {'': 4 * 'N'}
            for s in (''.join(i) for i in itertools.product(*(4 * ['ATGC']))):
                if rc(s) not in seq_map[4] or s not in seq_map[4]:
                    seq_map[4][s] = s
                    seq_map[4][rc(s)] = s
            seq_map[3] = {'': 3 * 'N'}
            for s in (''.join(i) for i in itertools.product(*(3 * ['ATGC']))):
                seq_map[3][s] = s
            seq_map[2] = {'': 2 * 'N'}
            for s in (''.join(i) for i in itertools.product(*(2 * ['ATGC']))):
                seq_map[2][s] = s

        #subsequences generator
        abun = {2: {'NN': 0}, 3: {'NNN': 0}, 4: {'NNNN': 0}}
        for l in [2, 3, 4]:
            abun[l] = dict([(s, 0) for s in seq_map[l].values()])
            for ss in self.slid_win(l):
                abun[l][seq_map[l].get(ss, l * 'N')] += 1
        zscore = {}
        for tet, f in abun[4].items():
            if f != 0 and tet != 'NNNN':
                n23 = abun[2][tet[1:3]]
                if n23 != 0:
                    n123 = abun[3][tet[:3]]
                    n234 = abun[3][tet[1:]]
                    e = n123 * n234 / n23
                    v = e * (n23 - n123) * (n23 - n234) / n23 ** 2
                else:
                    e, v = 0, 0
                n23i = abun[2][rc(tet[1:3])]
                if n23i != 0:
                    n123i = abun[3][rc(tet[:3])]
                    n234i = abun[3][rc(tet[1:])]
                    ei = n123i * n234i / n23i
                    vi = ei * (n23i - n123i) * (n23i - n234i) / n23i ** 2
                else:
                    ei, vi = 0, 0
                tv = sqrt(sqrt(v ** 2 + vi ** 2))
                if tv == 0:
                    tv = 1e-8
                zscore[tet] = (f - e - ei) / tv
            else:
                zscore[tet] = 0
        return zscore

    def mass(self):
        RNA_mass_table = {'A': 347.22, 'U': 324.18, 'G': 363.22, 'C': 323.20}
        DNA_mass_table = {'A': 331.22, 'T': 320.19, 'G': 347.22, 'C': 307.20}
