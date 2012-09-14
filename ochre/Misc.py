def letter_library(letter, alphabet):
    if alphabet == 'amino acid':
        lookup = {'R': 'Arginine',
                  'H': 'Histidine',
                  'K': 'Lysine',
                  'D': 'Aspartic Acid',
                  'E': 'Glutamic Acid',
                  'S': 'Serine',
                  'T': 'Threonine',
                  'N': 'Asparagine',
                  'Q': 'Glutamine',
                  'C': 'Cysteine',
                  'G': 'Glycine',
                  'P': 'Proline',
                  'A': 'Alanine',
                  'V': 'Valine',
                  'I': 'Isoleucine',
                  'L': 'Leucine',
                  'M': 'Methionine',
                  'F': 'Phenylalanine',
                  'Y': 'Tyrosine',
                  'W': 'Tryptophan',
                  'B': 'Asparagine or Aspartic Acid (N or D)',
                  'J': 'Isoleucine or Leucine (I or L)',
                  'O': 'Pyrrolysine',
                  'U': 'Selenocysteine',
                  'Z': 'Glutamine or Glutamic Acid (Q or E)',
                  'X': 'Any Residue',
                  '*': 'Stop',
                  '-': 'Gap'}
    elif alphabet == 'nucleic acid':
        lookup = {'A': 'Adenine',
                  'C': 'Cytosine',
                  'G': 'Guanine',
                  'T': 'Thymine',
                  'N': 'Any Residue',
                  'R': 'Purine (A or G)',
                  'Y': 'Pyrimide (C or T or U)',
                  'U': 'Uracil',
                  '-': 'Gap'}
    return lookup[letter]


def seq3_table(letter):
    lookup = {'R': 'Arg', 'H': 'His', 'K': 'Lys', 'D': 'Asp',
              'E': 'Glu', 'S': 'Ser', 'T': 'Thr', 'N': 'Asn',
              'Q': 'Gln', 'C': 'Cys', 'G': 'Gly', 'P': 'Pro',
              'A': 'Ala', 'V': 'Val', 'I': 'Ile', 'L': 'Leu',
              'M': 'Met', 'F': 'Phe', 'Y': 'Tyr', 'W': 'Trp',
              'B': 'Asx', 'J': 'Xle', 'O': 'Pyl', 'U': 'Sel',
              'Z': 'Glx', '*': 'Ter'}
    return lookup.get(letter, 'Xaa')


def tTables(tableName='standard'):
    #TODO: codes from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c
    if tableName == 'standard':
        return {'F': ['TTT','TTC'],
                'L': ['TTA','TTG','CTT','CTC','CTA','CTG'],
                'I': ['ATT','ATC','ATA'],
                'M': ['ATG'],
                'V': ['GTT','GTC','GTA','GTG'],
                'S': ['TCT','TCC','TCA','TCG','AGT','AGC'],
                'P': ['CCT','CCC','CCA','CCG'],
                'T': ['ACT','ACC','ACA','ACG'],
                'A': ['GCT','GCC','GCA','GCG'],
                'Y': ['TAT','TAC'],
                'H': ['CAT','CAC'],
                'Q': ['CAA','CAG'],
                'N': ['AAT','AAC'],
                'K': ['AAA','AAG'],
                'D': ['GAT','GAC'],
                'E': ['GAA','GAG'],
                'C': ['TGT','TGC'],
                'W': ['TGG'],
                'R': ['CGT','CGC','CGA','CGG','AGA','AGG'],
                'G': ['GGT','GGC','GGA','GGG'],
                '*': ['TAA','TAG','TGA']}
