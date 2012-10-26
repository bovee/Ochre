#!/usr/bin/env python
import sys
import os.path as op
sys.path.append(op.join(op.dirname(op.realpath(__file__)), '..'))
from ochre import seqlist


def summary(seqs, outfile):
    slen, bps = len(seqs), sum(len(s) for s in seqs)
    outfile.write('Sequences: ' + str(slen) + '\n')
    outfile.write('Basepairs: ' + str(bps) + '\n\n')
    outfile.write('Shortest: ' + str(min(len(s) for s in seqs)) + ' bps\n')
    outfile.write('N90: ' + str(seqs.n(0.9)) + ' bps\n')
    outfile.write('Average: ' + str(bps / slen) + ' bps\n')
    outfile.write('N50: ' + str(seqs.n(0.5)) + ' bps\n')
    outfile.write('N10: ' + str(seqs.n(0.1)) + ' bps\n')
    outfile.write('Longest: ' + str(max(len(s) for s in seqs)) + ' bps\n')


def bycontig(seqs, outfile, dl=',', kmer=100, megan=None, cov=None):
    if megan is None:
        outfile.write(dl.join(['length','gc','coverage'])+'\n')
    else:
        n2clade = dict(line.split(',') for line in megan)
        outfile.write(dl.join(['length','gc','coverage','clade'])+'\n')

    if hasattr(cov, 'read'):
        n2cov = dict(line.split(',') for line in cov)

    for s in seqs:
        nm = s.name.split(' ')[0]
        ls = len(s)
        if hasattr(cov, 'read'):
            cv = float(n2cov.get(nm, 0))
        elif cov == 'velvet':
            cv = float(s.name.split('_')[5]) * len(s) / int(s.name.split('_')[3])
        elif cov == 'idba':
            #this is an underestimate
            cv = float(s.name.split('_')[-1]) * ls / (ls - kmer + 1)
        if megan is None:
            outfile.write(str(ls) + dl + str(s.gc()) + dl + str(cv) + '\n')
        else:
            cld = n2clade.get(s.name.split(' ')[0], 'None').strip()
            outfile.write(str(ls) + dl + str(s.gc()) + dl + str(cv) + dl + cld + '\n')


def gc(seqs, outfile):
    outfile.write('gc\n')
    outfile.write('\n'.join(str(s.gc()) for s in seqs))


def tetra(seqs, outfile, dl=','):
    import itertools

    def rc(seq):  # reverse complement
        invert_table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(invert_table.get(i, 'N') for i in seq[::-1])

    seq_map = {'': 4 * 'N'}
    for s in (''.join(i) for i in itertools.product(*(4 * ['ATGC']))):
        if rc(s) not in seq_map or s not in seq_map:
            seq_map[s] = s
            seq_map[rc(s)] = s

    # write out the header
    srted_vals = list(set(seq_map.values()))
    srted_vals.sort()
    outfile.write(dl.join(['Gene'] + srted_vals) + '\n')

    for s in seqs:
        frq = s.nuc_freqs(seq_map)
        outfile.write(dl.join([s.name] + \
          [str(frq[i]) for i in srted_vals]))
        outfile.write('\n')
        outfile.flush()


def tetraz(seqs, outfile, dl=','):
    import itertools

    def rc(seq):  # reverse complement
        invert_table = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(invert_table.get(i, 'N') for i in seq[::-1])

    seq_map = {}
    seq_map[4] = {'': 4 * 'N'}
    for s in (''.join(i) for i in itertools.product(*(4 * ['ATGC']))):
        if rc(s) not in seq_map[4] or s not in seq_map[4]:
            seq_map[4][s] = s
            seq_map[4][rc(s)] = s
    seq_map[3] = {'': 3 * 'N'}
    for s in (''.join(i) for i in itertools.product(*(3 * ['ATGC']))):
        seq_map[3][s] = s
    seq_map[2] = {'': 2 * 'N'}
    for s in (''.join(i) for i in itertools.product(*(2 * ['ATGC']))):
        seq_map[2][s] = s

    # write out the header
    srted_vals = list(set(seq_map[4].values()))
    srted_vals.sort()
    outfile.write(dl.join(['Gene'] + srted_vals) + '\n')

    for s in seqs:
        frq = s.tetra_zscore(seq_map)
        outfile.write(dl.join([s.name] + \
          [str(frq[i]) for i in srted_vals]))
        outfile.write('\n')
        outfile.flush()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description= \
      'Calculate statistics for different types of sequence files.')
    parser.add_argument('type', choices=('summary', 'bycontig', 'gc', 'tetra', 'tetraz'), \
      nargs='?', default='summary', \
      help='Type of statistical analysis to run.')
    parser.add_argument('infile', \
      help='Name of the sequence file.')
    parser.add_argument('--outfile', '-o', type=argparse.FileType('w'), \
      default=sys.stdout, help='Name of the file to output.')
    parser.add_argument('--kmer', '-k', type=int, default=0, \
      help='For Velvet analysis, the k-mer length. For IDBA, the longest k-mer used.')
    parser.add_argument('--first', type=int, \
      help='How many sequences to analyse, starting with the first.')
    parser.add_argument('--megan', type=argparse.FileType('r'), default=None, \
      help='A CSV file from Megan with "read_name, taxon_name".')
    parser.add_argument('--coverage', type=argparse.FileType('r'), default=None, \
      help='A CSV file with "read_name, coverage".')
    parser.add_argument('--assembler', choices=('idba', 'velvet', 'none'), \
      default='none', help='The assembler used to produce the contigs.')
    #TODO: loose indexing for gc, bycontig and tetra?
    args = parser.parse_args()
    infile = open(args.infile, 'rb')
    seqs = seqlist(infile)
    if args.first is not None:
        seqs = seqs[:args.first]
    if args.type == 'bycontig':
        if args.coverage is not None:
            cov = args.coverage
        else:
            cov = args.assembler
        bycontig(seqs, args.outfile, kmer=args.kmer, megan=args.megan, cov=cov)
    elif args.type == 'summary':
        summary(seqs, args.outfile)
    elif args.type == 'gc':
        gc(seqs, args.outfile)
    elif args.type == 'tetra':
        tetra(seqs, args.outfile)
    elif args.type == 'tetraz':
        tetraz(seqs, args.outfile)
