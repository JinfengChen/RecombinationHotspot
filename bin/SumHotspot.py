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
python SumHotspot.py --input ../input/MPR.cross.uniq --bin ../input/MPR.geno.bin.uniq

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

'''Kick out chr from full chrlist if this line (rank) does not contain genetic distance of this chr'''
def listvalidchr(chrlist, binn, rank):
    copylist = chrlist[:]
    #print copylist
    for i in range(len(binn)):
        if int(rank) > int(binn[i]):
            #print i
            copylist.remove(int(i)+1)
    return copylist
   
 
'''Get genetic distance from map file'''
def Genetic_distance(mapfile):
    chrlist   = range(1,13)
    binn      = []
    data      = []
    gdistance = {}
    markercm  = defaultdict(lambda : defaultdict(float))
    markername = defaultdict(list) 
    with open(mapfile, 'r') as mapfh:
        s0   = re.compile(r'^-l')
        s1   = re.compile(r'^-Number')
        s2   = re.compile(r'^-b MarkerNames')
        s3   = re.compile(r'^-e MarkerNames')
        '''1 1 0100222046\n1 2 0100609676'''
        s4   = re.compile(r'^(\d+) (\d+) (\d+)')
        flag = 0
        for line in mapfh:
            if s1.search(line):
                line = line.rstrip()
                unit = re.split(r'\s+',line)
                binn = unit[2:]
                #print binn 
            if s0.search(line):
                data.append(line)
            if s2.search(line):
                flag = 1
            if s3.search(line):
                flag = 0
            if flag == 1 and s4.search(line):
                m = s4.search(line)
                chrname = int (m.groups(0)[0])
                marker  = m.groups(0)[2]
                marker  = int(marker[2:])
                markername[chrname].append(marker)

    for line in data:
        line = line.rstrip()
        unit = re.split(r'\s+',line)
        distanceline = unit[3:]
        '''skip the first line, which do not have genetic distance. For all genetic distance, we use marker1 to reference the physical position of the genetic distance. As for recombination bin, marker1 is the boundary of the first bin of the genetic interval, which close to recombination breakpoint'''
        if int(unit[1]) == int(0):
            continue
        #print unit
        newchrlist = listvalidchr(chrlist, binn, unit[1])
        #print newchrlist
        for i in range(len(distanceline)):
            gdistance.setdefault(newchrlist[i],[]).append(distanceline[i])

    for c in sorted(gdistance.keys()):
        #print c
        '''skip last marker, because the last one do not have genetic distance'''
        for i in range(len(gdistance[c])-1):
            #print markername[c][i],gdistance[c][i]
            markercm[int(c)][int(markername[c][i])]= float(gdistance[c][i])   
    return markercm

'''
Find recombination hotspot using number of recombinants in marker interval.
These significantly larger than genome average regions will be hotspot
'''
def hotspot_recombinant (binfile, markerdata):
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
        lastmkcm  = 0
        lastmkmid = 0
        for p in sorted (breakpoint[c].keys()):
            #print c, p, breakpoint[c][p]
            '''genetic distance start from the second marker and compared with the first marker of each chromosome'''
            if lastmkpos == 0:
                mid       = (int(p)-lastmkpos)/2
                lastmkpos = int(p)
                lastmkcm  = 0
                lastmkmid = mid
                continue
            else:
                '''physical distance between two bin is the distance from the midpoint of two bin'''
                mid       = float(lastmkpos) + (int(p)-lastmkpos)/2
                physicald = int(p) - int(lastmkpos)
                #physicald = float(mid) - float(lastmkmid)
                #print int(p), lastmkpos, mid, lastmkmid, physicald
                '''recombiantion frequency: r = R/2(1-R). Refrence A High-Resolution Map of Arabidopsis Recombinant Inbred Lines by Whole-Genome Exon Array Hybridization'''
                R         = float(breakpoint[c][p])/float(samplen)
                recfreq   = R/2.00*(1-R)
                '''Haldane mapping function: m = -(1/2) ln(1-2c), cM = 100 * m'''
                geneticdh  = -0.5 * float(np.log(1 - 2 * recfreq)) * 100
                '''Kosamibi mapping function: m = (1/4) ln((1+2c)/(1-2c))'''
                geneticdk = 0.25 * float(np.log((1+2*recfreq)/(1-2*recfreq))) * 100
                geneticd  = markerdata[int(c)][lastmkpos]
                markercm  = float(lastmkcm) + geneticd
                #lastmkpos = int(p)
                #lastmkcm  = markercm
                #print p, lastmkpos, physicald, breakpoint[c][p], samplen, geneticd, win
                geneticrate = geneticd/(float(physicald)/1000000.00)
                mean, pvalue = binominal_hotspot(breakpoint[c][p],physicald)
                print >> fh, int(c),int(lastmkpos),breakpoint[c][p],geneticd,physicald, geneticrate, mean, pvalue, markercm
                lastmkpos = int(p)
                lastmkcm  = markercm
                lastmkmid = mid
    fh.close()

    Rcmd ='''
data <- read.table("HEG4vsNB.Bimonimal.Hotspot.table")
chrlen <- c(43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856)
cent   <- matrix(c(16701176,17133774, 13579371,13757240, 19541173,19633542, 9757057,9880597, 12458897,12555251, 15424916,15490360, 11961237,12281215, 12920962,13840070, 2750852,2981329, 8100966,8178267, 12046080,12329552, 11772715,12070903),nrow=2,ncol=12)

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
    rect(cent[1,i],min(chrdata[,6]),cent[2,i],max(chrdata[,6]),border=NA,col=rgb(0.3,1,1,alpha=0.5))
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

pdf("HEG4vsNB.Genetic.Physical.pdf",12,18)
par(mfrow=c(6,2))
par(mar=c(4,4,3,3))
for (i in 1:12){
    chr = data[,1] == i
    chrdata <- data[chr,]
    color <- rep('black',length(chr))
    plot(chrdata[,2],chrdata[,9],xlim=c(0,45000000),col=color,xlab='',ylab='',axes=FALSE)
    rect(cent[1,i],min(chrdata[,9]),cent[2,i],max(chrdata[,9]),border=NA,col=rgb(0.3,1,1,alpha=0.5))
    ##x axis
    atx <- c(seq(0,ceiling(chrlen[i]/1000000),by=5),ceiling(chrlen[i]/1000000))
    axis(1,at=atx*1000000,labels=atx)
    mtext(paste("Chr",i," Physical Position (Mb)", sep=''), side=1, at = 15 * 1000000, line=2.5, cex = 0.8)
    ##y axis
    aty <- c(seq(0,max(chrdata[,9]),by= ceiling(max(chrdata[,9]-0)/5/10)*10),ceiling(max(chrdata[,9])))
    axis(2,at=aty,labels=aty)
    mtext("Genetic Position (cM)", side=2, at = ceiling(max(chrdata[,9]))/2, line=2.5, cex = 0.8)
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
                '''p1 is NB and p2 is HEG4'''
                p1 = genotype[c][p][0]
                p2 = genotype[c][p][1]
                pvalue = chisquare_test([p1,p2],[int(samplen/2),int(samplen/2)])
                print >> markerfh, int(c),int(p),p1,p2,float(p1)/float(p2),int(samplen/2),pvalue
    
    Rcmd ='''
data <- read.table("HEG4vsNB.Marker.Segregation.table")
chrlen <- c(43270923,35937250,36413819,35502694,29958434,31248787,29697621,28443022,23012720,23207287,29021106,27531856)
cent   <- matrix(c(16701176,17133774, 13579371,13757240, 19541173,19633542, 9757057,9880597, 12458897,12555251, 15424916,15490360, 11961237,12281215, 12920962,13840070, 2750852,2981329, 8100966,8178267, 12046080,12329552, 11772715,12070903),nrow=2,ncol=12)

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
    rect(cent[1,i],0.5,cent[2,i],2,border=NA,col=rgb(0.3,1,1,alpha=0.5))
    ##x axis
    atx <- c(seq(0,ceiling(chrlen[i]/1000000),by=5),ceiling(chrlen[i]/1000000))
    axis(1,at=atx*1000000,labels=atx)
    mtext(paste("Chr",i," Physical Position (Mb)", sep=''), side=1, at = 15 * 1000000, line=2.5, cex = 0.8)
    ##y axis
    aty <- seq(0.5,2,by=0.5)
    axis(2,at=aty,labels=aty)
    mtext("Segregation Ratio (NB/HEG4)", side=2, at = 1.2, line=2.5, cex = 0.8)
}
dev.off()
'''
    with open('HEG4vsNB.Marker.Segregation.R','w') as Rscript:
        Rscript.write(Rcmd)
    os.system('cat HEG4vsNB.Marker.Segregation.R | R --slave') 


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-b', '--bin')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    
    binfile = args.bin
    mapfile = args.input + '.map'
    markercm= Genetic_distance(mapfile)
    hotspot_recombinant(binfile, markercm)
    #marker_segregation(binfile) 
    #chisquare_hotspot(10,108) 
    #chisquare_hotspot(20,98) 
    #chisquare_hotspot(30,88)   

if __name__ == '__main__':
    main()

