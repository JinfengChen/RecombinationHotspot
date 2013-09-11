#!/opt/Python/2.7.3/bin/python
import sys
import re
import os
from collections import defaultdict
import numpy as np
from numpy import *
from scipy.stats.mstats import chisquare
from scipy.stats import fisher_exact
from scipy.stats import binom_test
import argparse

def usage():
    message='''
python SumHotspot.py --input ../input/MPR.geno.bin.uniq

    '''
    print message

'''
Method refered from Drouaud et al, 2006 Genome Research:
Variation in crossing-over rates across chromosome 4 of Arabidopsis thaliana reveals the presence of meiotic recombination hot spots

P is the probability of having a CO in a specific position of one chromosome
L is the length of interval between two marker
p = P*L is the probability of CO in interval of length L
k = rec is the recombinant (CO) that observed in the interval
V = sample is the number of chromosomes for the interval
binomial distribution: B(V, P*L)
'''
def binominal_hotspot(rec, length, sample=118, recombinant=3684, gsize=373245519):
    P = float(recombinant)/(float(sample) * float(gsize))
    L = length
    p = P * L
    k = rec
    V = sample
    #print k, V, P, L, p
    mean = int (p * sample)
    pvalue = binom_test(k, V, p) 
    return mean, pvalue

def chisquare_test(observed0, expected0):
    observed = np.array(observed0)
    expected = np.array(expected0)
    if min(observed0) > 5:
        a = chisquare(observed, expected)
    else:
        a = fisher_exact([observed, expected])
    return a[1]
    

'''
Find recombination hotspot using number of recombinants in marker interval.
These significantly larger than genome average regions will be hotspot
'''
def hotspot_recombinant (binfile):
    breakpoint = defaultdict(lambda : defaultdict(int))
    genotype   = defaultdict(lambda : defaultdict(str))
    samplen    = 0
    s = re.compile(r'^(\d{2})(\d+)')
    with open (binfile,'r') as binfh:
        '''Read header and store rils name as numberic for output'''
        rilsline = binfh.readline()
        rilsline = re.sub(r'"',r'',rilsline)
        rilsline = re.sub(r'GN',r'',rilsline)
        rilsline = rilsline.rstrip()
        rils     = rilsline.split('\t')
        samplen  = len(rils) - 1
        '''read matrix of recombination bin, find breakpoint for each RILs each Chromosome'''
        for line in binfh:
            line = re.sub(r'"',r'',line)
            m = s.search(line)
            if m:
                chro = m.groups(0)[0]
                pos  = m.groups(0)[1]
                line = line.rstrip()
                bins = line.split('\t')
                for i in range (1, len(bins)):
                    '''test if breakpoint is found, which mean 0 to 1 ot 1 to 0 transient'''
                    tag = 0
                    '''if previous genotype is defined, compare and let tag = 1 if breakpoint found. the default value for dict is '' '''
                    if genotype[rils[i]][chro] != '': 
                        tag = 0 if float (bins[i]) == float (genotype[rils[i]][chro]) else 1
                    '''update the genotype for the last record, which will used to compare in the next loop'''
                    genotype[rils[i]][chro] = bins[i]

                    '''update breakpoint number for data matirx'''
                    breakpoint[chro][pos] += tag
                    #print bins[0], rils[i], breakpoint[bins[0]]
    '''Output recombinant and non-recombinant for each windows'''
    recom_win_rec    = defaultdict(lambda : defaultdict(int))
    recom_win_nonrec = defaultdict(lambda : defaultdict(int))
    recom_win_mk = defaultdict(lambda : defaultdict(int))
    recom_win_recrate= defaultdict(lambda : defaultdict(list)) 
    recom_win_physicald= defaultdict(lambda : defaultdict(list))
    fh = open('HEG4vsNB.Bimonimal.Hotspot.table','w')
    for c in sorted (breakpoint.keys()):
        lastmkpos = 0
        for p in sorted (breakpoint[c].keys()):
            #print c, p, breakpoint[c][p]
            '''genetic distance start from the second marker and compared with the first marker of each chromosome'''
            if lastmkpos == 0:
                lastmkpos = int(p)
                continue
            else:
                physicald = int(p) - int(lastmkpos)
                geneticd  = 100.00 * breakpoint[c][p]/samplen
                lastmkpos = int(p)
                #print p, lastmkpos, physicald, breakpoint[c][p], samplen, geneticd, win
                geneticrate = geneticd/(float(physicald)/1000000.00)
                mean, pvalue = binominal_hotspot(breakpoint[c][p],physicald)
                print >> fh, int(c),int(p),breakpoint[c][p],geneticd,physicald, geneticrate, mean, pvalue
    fh.close()

    Rcmd ='''
data <- read.table("HEG4vsNB.Bimonimal.Hotspot.table")
chrlen <- c(43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856)
pdf("HEG4vsNB.Bimonimal.Hotspot.pdf",12,18)
par(mfrow=c(6,2))
par(mar=c(4,4,3,3))
for (i in 1:12){
    chr = data[,1] == i
    chrdata <- data[chr,]
    significant <- chrdata[,8] <= 0.05
    hotspot <- chrdata[,3] > chrdata[,7]
    hotspotsig <- intersect(which(significant),which(hotspot))
    color <- rep('black',length(chr))
    color[hotspotsig] <- 'red'
    plot(chrdata[,2],chrdata[,6],xlim=c(0,45000000),col=color,xlab='',ylab='',axes=FALSE)
    ##x axis
    atx <- c(seq(0,ceiling(chrlen[i]/1000000),by=5),ceiling(chrlen[i]/1000000))
    axis(1,at=atx*1000000,labels=atx)
    mtext(paste("Chr",i," Physical Position (Mb)", sep=''), side=1, at = 15 * 1000000, line=2.5, cex = 0.8)
    ##y axis
    aty <- c(seq(0,max(chrdata[,6]),by= ceiling(max(chrdata[,6]-0)/5/10)*10),ceiling(max(chrdata[,6])))
    axis(2,at=aty,labels=aty)
    mtext("Recombination Rate (cM/Mb)", side=2, at = ceiling(max(chrdata[,6]))/2, line=2.5, cex = 0.8)
}
dev.off()
'''
    with open('HEG4vsNB.Binominal.Hotspot.R','w') as Rscript:
        Rscript.write(Rcmd)
    os.system('cat HEG4vsNB.Binominal.Hotspot.R | R --slave') 


'''
Marker segregation
'''
def marker_segregation(binfile):
    genotype   = defaultdict(lambda : defaultdict(lambda : defaultdict(int)))
    samplen    = 0
    s = re.compile(r'^(\d{2})(\d+)')
    with open (binfile,'r') as binfh:
        '''Read header and store rils name as numberic for output'''
        rilsline = binfh.readline()
        rilsline = re.sub(r'"',r'',rilsline)
        rilsline = re.sub(r'GN',r'',rilsline)
        rilsline = rilsline.rstrip()
        rils     = rilsline.split('\t')
        samplen  = len(rils) - 1
        '''read matrix of recombination bin, find breakpoint for each RILs each Chromosome'''
        for line in binfh:
            line = re.sub(r'"',r'',line)
            m = s.search(line)
            if m:
                chro = m.groups(0)[0]
                pos  = m.groups(0)[1]
                line = line.rstrip()
                bins = line.split('\t')
                for i in range (1, len(bins)):
                    if int(bins[i]) == 1:
                        genotype[chro][pos][1] += 1
                    elif int(bins[i]) == 0:
                        genotype[chro][pos][0] += 1
                    #print chro, pos, bins[i],genotype[chro][pos][1],genotype[chro][pos][0]
                
    with open ('HEG4vsNB.Marker.Segregation.table','w') as markerfh:
        for c in sorted (genotype.keys()):
            for p in sorted (genotype[c].keys()):
                p1 = genotype[c][p][0]
                p2 = genotype[c][p][1]
                pvalue = chisquare_test([p1,p2],[int(samplen/2),int(samplen/2)])
                print >> markerfh, int(c),int(p),p1,p2,float(p1)/float(p2),int(samplen/2),pvalue
    
    Rcmd ='''
data <- read.table("HEG4vsNB.Marker.Segregation.table")
chrlen <- c(43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856)
pdf("HEG4vsNB.Marker.Segregation.pdf",12,18)
par(mfrow=c(6,2))
par(mar=c(4,4,3,3))
for (i in 1:12){
    chr = data[,1] == i
    chrdata <- data[chr,]
    significant <- chrdata[,7] <= 0.05
    markersig <- which(significant)
    color <- rep('black',length(chr))
    color[markersig] <- 'red'
    plot(chrdata[,2],chrdata[,5],xlim=c(0,45000000),ylim=c(0.5,2),col=color,xlab='',ylab='',axes=FALSE)
    ##x axis
    atx <- c(seq(0,ceiling(chrlen[i]/1000000),by=5),ceiling(chrlen[i]/1000000))
    axis(1,at=atx*1000000,labels=atx)
    mtext(paste("Chr",i," Physical Position (Mb)", sep=''), side=1, at = 15 * 1000000, line=2.5, cex = 0.8)
    ##y axis
    aty <- seq(0.5,2,by=0.5)
    axis(2,at=aty,labels=aty)
    mtext("Segregation Ratio", side=2, at = 1.2, line=2.5, cex = 0.8)
}
dev.off()
'''
    with open('HEG4vsNB.Marker.Segregation.R','w') as Rscript:
        Rscript.write(Rcmd)
    os.system('cat HEG4vsNB.Marker.Segregation.R | R --slave') 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    binfile = args.input
    #hotspot_recombinant(binfile)
    marker_segregation(binfile) 
    #chisquare_hotspot(10,108) 
    #chisquare_hotspot(20,98) 
    #chisquare_hotspot(30,88)   

if __name__ == '__main__':
    main()

