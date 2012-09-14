#Trimmomatic(, leading=3,trailing=3,slidingwin=(4,15),minlen=36)
from ochre.FileSeqList import PairedFileSeqList
from ochre.SeqList import SeqList
from ochre.External.External import app, run


def Trimmomatic(seqs, headcrop=None, leading=None, trailing=None, slidingwin=None, minlen=None, illuminaclip=None):
    if isinstance(seqs, PairedFileSeqList):
        tend = 'PE'
        f1, f2 = seqs._f1.get_file('fq'), seqs._f2.get_file('fq')
        flist = [f1, f2, \
                f1 + '.trim.prd.fq', f1 + '.trim.uprd.fq', \
                f2 + '.trim.prd.fq', f2 + '.trim.uprd.fq']
    elif isinstance(seqs, SeqList):
        tend = 'SE'
        f1 = seqs.get_file('fq')
        flist = [f1, f1 + '.trim.fq']
    else:
        raise TypeError('Need a PairedFileSeqList, FileSeqList or SeqList.')

    qual = '-phred33'  # TODO: determine this from the data
    c = [app('JAVA', 'java'), '-classpath', \
      app('TRIMMOMATIC', 'trimmomatic-0.22.jar'), \
      'org.usadellab.trimmomatic.Trimmomatic' + tend, qual]
    c += flist

    if illuminaclip is not None:
        c += ['ILLUMINACLIP:' + str(illuminaclip.get_file('fa')) + \
          ':2:40:15']  # sane defaults?
    if leading is not None:
        c += ['LEADING:' + str(leading)]
    if trailing is not None:
        c += ['TRAILING:' + str(trailing)]
    if slidingwin is not None:
        c += ['SLIDINGWINDOW:' + str(slidingwin[0]) + ':' + \
          str(slidingwin[1])]
    if headcrop is not None:
        c += ['HEADCROP:' + str(headcrop)]
    if minlen is not None:
        c += ['MINLEN:' + str(minlen)]
    run([c])
    if tend == 'PE':
        return PairedFileSeqList(f1 + '.trim.prd.fq',\
                                 f2 + '.trim.prd.fq')
    else:
        return SeqList(f1 + '.trim.fq')
