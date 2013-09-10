#!/opt/Python/2.7.3/bin/python
import sys
import re
import os
from collections import defaultdict
import numpy as np
from numpy import *
from scipy.stats.mstats import chisquare
from scipy.stats import fisher_exact
import argparse

def usage():
    message='''
python SumHotspot.py --input ../input/MPR.geno.bin.uniq

    '''
    print message

def chisquare_hotspot(rec, norec, sample=118, winsize=50000, binn=1614, recombinant=3684, gsize=372000000):
    meanbps_win = float(binn)/(gsize/winsize)
    totalrec_win= float(meanbps_win * float(sample))
    factor = 1
    meanrec_win = float(recombinant)/(gsize/winsize)
    '''if mean recombinant in win smaller than 1 time, we make it 1 and multiple total recombinant by factor that adapted to mean'''
    if meanrec_win < 1:
        factor = 1//meanrec_win
        meanrec_win  = 1
        totalrec_win = totalrec_win * factor
    
    meannorec_win = totalrec_win - meanrec_win
    observed = np.array([rec,norec])
    expected = np.array([meanrec_win,meannorec_win])
    a = chisquare(observed, expected)
    b = fisher_exact([observed, expected])
    return b[1]

'''
Find recombination hotspot using number of recombinants in 50 kb windows.
We have a genome-wide mean number of recombiants in 50 kb windows. So for every windows we do fisher exactly test.
These significantly larger than genome average regions will be hotspot
'''
def hotspot_recombinant (binfile,win=50000):
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
                index = int (p) // int (win)
                recom_win_physicald[c][index].append(physicald)
                recom_win_recrate[c][index].append(geneticrate)
                recom_win_rec[c][index] += breakpoint[c][p]
                recom_win_mk[c][index]  += 1
                recom_win_nonrec[c][index] = recom_win_mk[c][index] * samplen - recom_win_rec[c][index]
    for c in sorted (recom_win_rec.keys()):
        for i in sorted (recom_win_rec[c].keys()):
            #print c, i, recom_win_rec[c][i], recom_win_mk[c][i], recom_win_nonrec[c][i]
            pvalue = chisquare_hotspot(recom_win_rec[c][i], recom_win_nonrec[c][i])
            rate   = mean(recom_win_recrate[c][i])
            #print int(i) * int(win), pvalue
            if pvalue < 0.5:
                print c, int (i) * int(win), recom_win_rec[c][i], recom_win_mk[c][i], recom_win_nonrec[c][i], recom_win_recrate[c][i],recom_win_physicald[c][i], pvalue

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
    hotspot_recombinant(binfile) 
    chisquare_hotspot(10,108) 
    chisquare_hotspot(20,98) 
    chisquare_hotspot(30,88)   

if __name__ == '__main__':
    main()

