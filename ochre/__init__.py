"""
Ochre Sequence Handling Module
"""


def _add_to_path():
    import os
    import sys
    ochre._path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    sys.path.append(ochre._path)

#import things that are imported inside modules (to avoid circular conflicts)
#import ochre.FileSeqList.FileSeqList, ochre.FileSeqList.PairedFileSeqList
#import ochre.NucSeq.NASeq
#import ochre.PepSeq.PepSeq
#import ochre.FileFormats
from ochre.FileSeqList import FileSeqList, PairedFileSeqList
from ochre.NucSeq import NASeq
from ochre.PepSeq import PepSeq
import ochre.FileFormats.FASTA
import ochre.FileFormats.FASTQ
import ochre.FileFormats.Biopython

#for convenience import commonly used stuff
from ochre.Sequence import Seq as seq
from ochre.SeqList import SeqList as seqlist

_add_to_path()
