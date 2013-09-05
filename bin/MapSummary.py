#!/opt/Python/2.7.3/bin/python
import sys
sys.path.append("/rhome/cjinfeng/software/ProgramPython/module/lib")
from dictionary import dict_add
from numpy  import *
import re
import os
import argparse

def usage():
    message='''
python MapSummary.py --input ../input/MPR.cross.uniq --bin ../input/MPR.geno.bin.uniq  --output HEG4vsNB
The script read *.map and *.cro as input, parse and output general summary of the map.
--input: input is a prefix of .cro and .map. ../input/MPR.cross.uniq will stand for ../input/MPR.cross.uniq.cro and ../input/MPR.cross.uniq.map
--bin:   marker-RILs matrix of uniq bin map
--output: output is then prefix of output file. Optional, Defaule is HEG4vsNB.

    '''
    print message

'''Kick out chr from full chrlist if this line (rank) does not contain genetic distance of this chr'''
def listvalidchr(chrlist, binn, rank):
    copylist = chrlist[:]
    #print copylist
    for i in range(len(binn)):
        if int(rank) > int(binn[i]):
            #print i
            copylist.remove(int(i)+1)
    return copylist

'''Summary recombination bin, genetic distance, recombiantion rate for each chromosome'''
def Chr_breakpoint(mapfile,chrbp,recombinant):
    chrlist   = range(1,13)
    binn      = []
    data      = []
    gdistance = {}
    chrlen    = {1:43270923, 2:35937250, 3:36413819, 4:35502694, 5:29958434, 6:31248787, 7:29697621, 8:28443022, 9:23012720, 10:23207287, 11:29021106, 12:27531856} 
    with open(mapfile, 'r') as mapfh:
        s0   = re.compile(r'^-l')
        s1   = re.compile(r'^-Number')
        for line in mapfh:
            if s1.search(line):
                line = line.rstrip()
                unit = re.split(r'\s+',line)
                binn = unit[2:]
                print binn 
            if s0.search(line):
                data.append(line)
    for line in data:
        line = line.rstrip()
        unit = re.split(r'\s+',line)
        distanceline = unit[3:]
        if int(unit[1]) == int(0):
            continue
        #print unit 
        newchrlist = listvalidchr(chrlist, binn, unit[1])
        #print newchrlist
        for i in range(len(distanceline)):
            gdistance.setdefault(newchrlist[i],[]).append(distanceline[i])
    #print "\n".join(gdistance[1])
    chrfh = open (chrbp,'w')
    chrfh.write('Chr\t#Bin\t#Recombinant\tMean genetic distance (cM)\tTotal genetic distance (cM)\tPhysical length (Mb)\tRecombination Rate (cM/Mb)\n')
    sumall = {}
    distance = []
    for c in sorted (gdistance.keys()):
        distance.extend(gdistance[c])
        gdistance[c] = map (float,gdistance[c])
        cmean   = "%.2f" %mean(gdistance[c])
        cmedian = "%.2f" %median(gdistance[c])
        ctotal  = sum(gdistance[c])
        ctotal  = "%.2f" %ctotal
        cmin    = min(gdistance[c])
        cmax    = max(gdistance[c])
        clen    = "%.2f" %(float(chrlen[c])/1000000.00)
        crate   = float(ctotal)/float(clen)
        crate   = "%.2f" %crate
        sumall['binn']   = sumall['binn'] + int(binn[int(c)-1]) if sumall.has_key('binn') else int(binn[int(c)-1])
        sumall['total']  = sumall['total'] + float(ctotal) if sumall.has_key('total') else float(ctotal)
        sumall['genome'] = sumall['genome'] + float(clen) if sumall.has_key('genome') else float(clen)
        print c, cmean, cmedian, cmin, cmax, ctotal 
        crecom  = str(recombinant[int(c)-1])
        chrfh.write('Chr'+str(c)+'\t'+binn[int(c)-1]+'\t'+crecom+'\t'+str(cmean)+'\t'+str(ctotal)+'\t'+str(clen)+'\t'+str(crate)+'\n')
    sumall['rate'] = float(sumall['total'])/float(sumall['genome'])
    sumall['rate'] = "%.2f" %sumall['rate']
    sumall['genome'] = "%.2f" %sumall['genome']
    sumall['total'] = "%.2f" %sumall['total']
    sumall['distance'] = "%.2f" %mean(map(float, distance))
    chrfh.write('Total'+'\t'+str(sumall['binn'])+'\t'+str(recombinant[12])+'\t'+sumall['distance']+'\t'+sumall['total']+'\t'+sumall['genome']+'\t'+sumall['rate']+'\n')
    chrfh.close()

'''Summary recombination bin for each RIL on each chromosome'''
def RIL_breakpoint(binfile,rilsbp):
    genotype={}
    breakpoint={}
    s = re.compile(r'^(\d{2})')
    with open (binfile,'r') as binfh:
        '''Read header and store rils name as numberic for output'''
        rilsline = binfh.readline()
        rilsline = re.sub(r'"',r'',rilsline)
        rilsline = re.sub(r'GN',r'',rilsline)
        rilsline = rilsline.rstrip()
        rils     = rilsline.split('\t')
        '''read matrix of recombination bin, find breakpoint for each RILs each Chromosome'''
        for line in binfh:
            line = re.sub(r'"',r'',line)
            m = s.search(line)
            if m:
                chro = m.groups(0)[0]
                line = line.rstrip()
                bins = line.split('\t')
                for i in range (1, len(bins)):
                    '''test if breakpoint is found, which mean 0 to 1 ot 1 to 0 transient'''
                    tag = 0 
                    if genotype.has_key(rils[i]):
                        chrgenotype = genotype[rils[i]]
                        if chrgenotype.has_key(chro):
                            tag = 0 if float (bins[i]) == float (genotype[rils[i]][chro]) else 1
                    genotype   = dict_add(genotype, rils[i], chro, bins[i])

                    '''update breakpoint number for data matirx'''
                    newbp = 0
                    if breakpoint.has_key(rils[i]):
                        chrbreakpoint = breakpoint[rils[i]]
                        if chrbreakpoint.has_key(chro):
                            newbp = breakpoint[rils[i]][chro] + tag        
                    breakpoint = dict_add(breakpoint, rils[i], chro, newbp)

                    #print(rils[i], chro) 
                    #print("Tag", tag)
                    #print("Bin", bins[i])
                    #print("newbp", newbp)
                    #print("genotype", genotype[rils[i]][chro])
                    #print("breakpoint", breakpoint[rils[i]][chro])
    '''Output chromosome matrix of breakpoint number for each RIL'''
    rilsbpfh = open (rilsbp,'w')
    rilsheader = ''
    rilsmatrix = ''
    chrsum = {}
    for r in sorted (map (int, breakpoint.keys())):
        r = str(r)
        total = 0
        rilsbpdata=[]
        rilshead  =['RILs']
        rilsbpdata.append('GN' + r)
        breakpointc = breakpoint[r]
        for c in sorted (breakpointc.keys()):
            chromosome = int (c)
            chrsum[chromosome] = chrsum[chromosome] + int(breakpointc[c]) if chrsum.has_key(chromosome) else int(breakpointc[c])
            total += breakpointc[c]
            rilsbpdata.append(breakpointc[c])
            rilshead.append('Chr' + str(chromosome))
        rilsbpdata.append(total)
        rilshead.append('Total')
        rilsbpline = "\t".join(map(str,rilsbpdata))
        rilsheader = "\t".join(map(str,rilshead))
        rilsmatrix = rilsmatrix + rilsbpline + '\n'
     
    rilsbpfh.write(rilsheader + '\n')
    rilsbpfh.write(rilsmatrix)
    '''Sum for each chromosome and output total for each chromosome'''
    totalsum = []
    for c in sorted (chrsum.keys()):
        #print chrsum[c]
        totalsum.append(chrsum[c])
    totalsum.append(sum(totalsum))
    lastline = "\t".join(map (str, totalsum))
    rilsbpfh.write('Total\t'+str(lastline)+'\n')
    rilsbpfh.close()
    return totalsum
     
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
    if (args.output is None):
        args.output = 'HEG4svNB'       

    print 'INPUT file is' + args.input
    print 'OUTPUT file is' + args.output
    crossfile = args.input + '.cro'
    mapfile   = args.input + '.map'
    binfile   = args.bin
    chrbp     = args.output + '.Chr.Breakpoint.table'
    rilbp     = args.output + '.RIL.Breakpoint.table'
    summary   = args.output + '.General.Summary.table'
    
    recombinant = RIL_breakpoint(binfile,rilbp) 
    Chr_breakpoint(mapfile,chrbp,recombinant)
    #GeneralSum()
     
'''     
    for chr in range(1,13):
        datafile1  = args.output + '.Chr' + str(chr) +'.Datafile'
        markerline, markernum = datafilehead(mapfile,chr)
        datafilehaplotype(crossfile,datafile1,markerline,markernum,chr) 
'''
if __name__ == '__main__':
    main()

