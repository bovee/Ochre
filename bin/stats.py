#!/usr/bin/env python
import sys
import os.path as op
sys.path.append(op.join(op.dirname(op.realpath(__file__)),'..'))
from ochre import seqlist


def velvet(infile, outfile, kmer=1):
    #fh = open(outfile, 'w')
    seqs = seqlist(infile)
    print('length', 'gc', 'coverage', sep=',', file=outfile)
    data = []
    for s in seqs:
        cv = float(s.name.split('_')[5]) * len(s) / int(s.name.split('_')[3])
        print(len(s), s.gc(), str(cv), sep=',', file=outfile)

def normal(filename):
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description= \
            'Calculate statistics for different types of sequence files.')
    parser.add_argument('type', choices=('velvet','normal'), \
      nargs='?', default='normal', \
      help='Type of statistical analysis to run.')
    parser.add_argument('infile', \
      help='Name of the sequence file.')
    parser.add_argument('--outfile','-o', type=argparse.FileType('w'), \
      default=sys.stdout, \
      help='Name of the file to output.')
    parser.add_argument('--kmer','-k',type=int, default=0, \
      help='For Velvet analysis, the k-mer length.')
    args = parser.parse_args()
    infile = open(args.infile,'rb')
    if args.type == 'velvet':
        velvet(infile, args.outfile, args.kmer)
