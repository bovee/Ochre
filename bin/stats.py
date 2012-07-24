#!/usr/bin/env python
import sys
import os.path as op
sys.path.append(op.join(op.dirname(op.realpath(__file__)), '..'))
from ochre import seqlist


def velvet(infile, outfile, kmer=1):
    #fh = open(outfile, 'w')
    seqs = seqlist(infile)
    outfile.write('length,gc,coverage\n')
    for s in seqs:
        cv = float(s.name.split('_')[5]) * len(s) / int(s.name.split('_')[3])
        outfile.write(str(len(s)) + ',' + str(s.gc()) + ',' + str(cv) + '\n')


def idba(infile, outfile, kmer=100, megan=None):
    if megan is None:
        seqs = seqlist(infile)
        outfile.write('length,gc,coverage\n')
        for s in seqs:
            ls = len(s)
            cv = float(s.name.split('_')[-1]) * ls / (ls - kmer + 1)
            outfile.write(str(ls) + ',' + str(s.gc()) + ',' + str(cv) + '\n')
    else:
        seqs = seqlist(infile)
        name2clade = dict(line.split(',') for line in megan)
        outfile.write('length,gc,coverage,clade\n')
        for s in seqs:
            ls = len(s)
            cv = float(s.name.split('_')[-1]) * ls / (ls - kmer + 1)
            cld = name2clade(s.name.split(' ')[0],'None').strip()
            outfile.write(str(ls) + ',' + str(s.gc()) + ',' + str(cv) + ',' + cld + '\n')


def normal(infile, outfile):
    seqs = seqlist(infile)
    slen, bps = len(seqs), sum(len(s) for s in seqs)
    outfile.write('Sequences: ' + str(slen) + '\n')
    outfile.write('Basepairs: ' + str(bps) + '\n\n')
    outfile.write('Shortest: ' + str(min(len(s) for s in seqs)) + ' bps\n')
    outfile.write('N90: ' + str(seqs.n(0.9)) + ' bps\n')
    outfile.write('Average: ' + str(bps / slen) + ' bps\n')
    outfile.write('N50: ' + str(seqs.n(0.5)) + ' bps\n')
    outfile.write('N10: ' + str(seqs.n(0.1)) + ' bps\n')
    outfile.write('Longest: ' + str(max(len(s) for s in seqs)) + ' bps\n')


def gc(infile, outfile, n=None):
    seqs = seqlist(infile)
    outfile.write('gc\n')
    if n is None:
        outfile.write('\n'.join(str(s.gc()) for s in seqs))
    else:
        outfile.write('\n'.join(str(s.gc()) for s in seqs[:n]))


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description= \
      'Calculate statistics for different types of sequence files.')
    parser.add_argument('type', choices=('velvet', 'normal', 'idba', 'gc'), \
      nargs='?', default='normal', \
      help='Type of statistical analysis to run.')
    parser.add_argument('infile', \
      help='Name of the sequence file.')
    parser.add_argument('--outfile', '-o', type=argparse.FileType('w'), \
      default=sys.stdout, \
      help='Name of the file to output.')
    parser.add_argument('--kmer', '-k', type=int, default=0, \
      help='For Velvet analysis, the k-mer length. For IDBA, the sequence length.')
    parser.add_argument('--first', type=int, \
      help='How many sequences to analyse, starting with the first.')
    parser.add_argument('--megan', type=argparse.FileType('r'), default=0, \
      help='A CSV file from Megan with "read_name, taxon_name".')
    args = parser.parse_args()
    infile = open(args.infile, 'rb')
    if args.type == 'velvet':
        velvet(infile, args.outfile, args.kmer)
    elif args.type == 'normal':
        normal(infile, args.outfile)
    elif args.type == 'idba':
        idba(infile, args.outfile, args.kmer, args.megan)
    elif args.type == 'gc':
        gc(infile, args.outfile, args.first)
