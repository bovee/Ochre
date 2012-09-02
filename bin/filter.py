#!/usr/bin/env python
import sys
import os.path as op
sys.path.append(op.join(op.dirname(op.realpath(__file__)), '..'))
from ochre import seqlist


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description= \
      'Filter data from a sequence file into a new sequence file.')
    parser.add_argument('infile', \
      help='Name of the sequence file.')
    parser.add_argument('--outfile', '-o', type=argparse.FileType('w'), \
      default=sys.stdout, help='File to output.')
    parser.add_argument('--first', type=int, \
      help='How many sequences to analyse, starting with the first.')
    parser.add_argument('--maxlength', type=int, \
      help='Maximum sequence length.')
    parser.add_argument('--minlength', type=int, \
      default=0, help='Minimum sequence length.')
    parser.add_argument('--coverage', type=argparse.FileType('r'), \
      default=None, help='A CSV file with "read_name, coverage".')
    parser.add_argument('--maxcov', type=float, \
      help='Maximum coverage for each sequence.')
    parser.add_argument('--mincov', type=float, \
      default=0, help='Minimum coverage for each sequence.')
    args = parser.parse_args()
    infile = open(args.infile, 'rb')
    seqs = seqlist(infile)
    if args.first is not None:
        seqs = seqs[:args.first]

    if args.maxlength is not None:
        lseqs = [s for s in seqs if len(s) <= args.maxlength \
          and len(s) >= args.minlength]
    else:
        lseqs = [s for s in seqs if len(s) >= args.minlength]

    if args.coverage is not None:
        n2cov = dict(line.split(',') for line in args.coverage)
        if args.maxcov is not None:
            cseqs = [s for s in lseqs if n2cov.get(s.name.split(' ')[0], 0) >= args.mincov \
              and n2cov.get(s.name.split(' ')[0], 0) <= args.maxcov]
        else:
            cseqs = [s for s in lseqs if n2cov.get(s.name.split(' ')[0], 0) >= args.mincov]
    else:
        cseqs = lseqs

    if len(cseqs) > 0:
        seqlist(cseqs).write(args.outfile)
