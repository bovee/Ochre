"""
Ochre Sequence Handling Module
"""

def _add_to_path():
    import os, sys
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

#import things that are imported inside modules (to avoid circular conflicts)
from ochre.FileSeqList import FileSeqList, PairedFileSeqList
from ochre.NucSeq import NASeq
from ochre.PepSeq import PepSeq

#for convenience import commonly used stuff
from ochre.Sequence import Seq as seq
from ochre.SeqList import SeqList as seqlist

_add_to_path()
