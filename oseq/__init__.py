"""
Ochre Sequence Handling Module
"""

import os, sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__),'..')))

#import things that are imported inside modules (to avoid circular conflicts)
from oseq.FileSeqList import FileSeqList
from oseq.NucSeq import NASeq
from oseq.PepSeq import PepSeq
