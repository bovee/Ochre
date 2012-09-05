import ochre

with open('gc.csv', 'w') as f:
    for seq in ochre.seqlist('test.fa'):
        f.write(seq.name + ',' + str(seq.gc()))

# VERSUS

import Bio.SeqIO
import Bio.SeqUtils

with open('gc.csv', 'w') as f:
    for seq in Bio.SeqIO.parse('test.fa', 'fasta'):
        f.write(seq.name + ',' + str(Bio.SeqUtils.GC(seq)))
