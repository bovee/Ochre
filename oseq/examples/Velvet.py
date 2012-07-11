from oseq import seqlist


def get_velvet_stats(filename, kmer):
    seqs = seqlist(filename)
    print('Length', 'GC %', 'Velvet Length', 'Coverage')
    for s in seqs:
        ln = int(s.name.split('_')[3]) + kmer - 1
        cv = float(s.name.split('_')[5]) * ln / (ln - kmer + 1)
        print(len(s), s.gc(), str(ln), str(cv))
