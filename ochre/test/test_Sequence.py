import unittest
from oseq.Sequence import Seq


class Test(unittest.TestCase):
    """ Unit tests for the Sequence class. """

    def test_reverse(self):
        seq = 'AGTCAACT'
        rseq = Seq(seq).reverse().seq
        self.assertEqual(rseq, 'TCAACTGA')

if __name__ == '__main__':
    unittest.main()
