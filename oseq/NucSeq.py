from oseq.Sequence import Seq

class NASeq(Seq):
    def gc(self):
        #TODO: fails if R or Y present?
        return sum([1 for i in self.seq if i in 'GC'])/len(self.seq)

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
        invert_table = {'A':'T','U':'A','T':'A','C':'G','G':'C',
                        'R':'R', 'Y':'Y','N':'N','-':'-'}
        if is_rna:
            invert_table['A'] = 'U'
        return ''.join(invert_table[i] for i in seq)

    def melting_temp(self, method='basic'):
        #http://bioinformatics.oxfordjournals.org/content/21/6/711.long
        import collections
        bps = collections.Counter(self.seq.upper())
        if method == 'basic':
            return 64.9+41.0*(bps['G']+bps['C']-16.4) / \
              (bps['A'] + bps['T'] + bps['G'] + bps['C'])

    def find_orfs(self):
        raise NotImplementedError

    def intron_exon_scan(self):
        raise NotImplementedError

    def find_restrict_sites(self):
        raise NotImplementedError
    
    def find_repeats(self):
        raise NotImplementedError

    def design_primers(self):
        raise NotImplementedError

    def predict_struct(self):
        raise NotImplementedError
