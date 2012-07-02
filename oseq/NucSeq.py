from oseq.Sequence import Seq

class NASeq(Seq):
    def gc(self):
        #TODO: fails if R or Y present?
        return sum([1 for i in self.seq if i in 'GC'])/len(self.seq)

    def reverse(self):
        pass

    def invert(self):
        return self.__invert__()

    def _invert(self):
        invert_table = {'A':'T','U':'A','T':'A','C':'G','G':'C',
                        'R':'R', 'Y':'Y','N':'N','-':'-'}
        if self._is_rna():
            invert_table['A'] = 'U'
        return Seq(''.join(invert_table[i] for i in self.seq))

    def melting_temp(self):
        #http://en.wikipedia.org/wiki/DNA_melting
        melt_g =     {'AA':-4.26, 'TT':-4.26, 'AT':-3.67,
                      'TA':-2.50, 'CA':-6.12, 'GT':-6.09,
                      'CT':-5.40, 'GA':-5.51, 'CG':-9.07,
                      'GC':-9.36, 'GG':-7.66, 'CC':-7.66}
        raise NotImplementedError

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
        
class PairedNASeq(NASeq):
    def __init__(seq1,seq2):
        self.seq1 = seq1
        self.seq2 = seq2
