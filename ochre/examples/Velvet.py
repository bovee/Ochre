from ochre import seqlist
import numpy as np
import matplotlib.pyplot as plt


def get_velvet_stats(filename, kmer):
    fh = open('temp.csv', 'w')
    seqs = seqlist(filename)
    print('length', 'gc', 'coverage', sep=',', file=fh)
    data = []
    for s in seqs:
        cv = float(s.name.split('_')[5]) * len(s) / int(s.name.split('_')[3])
        print(len(s), s.gc(), str(cv), sep=',', file=fh)
    #library(ggplot2)
    #seqs <- read.csv('temp.csv')
    #attach(seqs)
    #pdf('graph.pdf')
    #ggplot(seqs, aes(x=gc,y=coverage,colour=length))+geom_point(pch=19,cex=0.5)+scale_color_gradient(low="red",high="blue",guide="colorbar",trans="log")+scale_y_log10()
    #dev.off()

def get_velvet_stats_matplotlib(filename, kmer):
    #fh = open('temp.csv', 'w')
    seqs = seqlist(filename)
    #print('Length', 'GC %', 'Coverage', sep=',', file=fh)
    data = []
    for s in seqs:
        cv = float(s.name.split('_')[5]) * len(s) / int(s.name.split('_')[3])
        #print(len(s), s.gc(), str(cv), sep=',', file=fh)
        data += [[len(s), s.gc(), cv]]
    data = np.array(data)
    plt.plot(data[:,2],data[:,1],c=data[:,0])
    plt.colorbar()

def stats_for_paper():
    #figure 1 - gc histogram (assemb. + raw)
    #figure 2 - cov (y) vs. contig length (x)
    #figure 3 - cov (y) vs. gc (all reads)
    #figure 4 - cov (y) vs. gc (reads > 3000 bp, colored by length)
    pass
