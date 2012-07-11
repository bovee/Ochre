from oseq.Sequence import Seq


class PepSeq(Seq):
    def pi(self):
        raise NotImplementedError

    def hydrophobicity(self):
        raise NotImplementedError
