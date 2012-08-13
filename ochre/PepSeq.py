from ochre.Sequence import Seq


class PepSeq(Seq):
    def to_nuc(self, table='standard', rna=False):
        from ochre.Misc import tTables
        tab = tTables(table)

        try:
            rseq = ''.join([tab[i][0] for i in self.seq])
        except KeyError:
            raise KeyError('Amino acid not found in table.')

        if rna:
            return Seq(rseq.replace('T', 'U'), seq_type='RNA')
        else:
            return Seq(rseq, seq_type='DNA')

    def to_pep(self, table='standard'):
        return self

    def pi(self):
        raise NotImplementedError

    def hydrophobicity(self):
        raise NotImplementedError
