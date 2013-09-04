#!/opt/Python/2.7.3/bin/python
import sys
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

def datafilehaplotype(cro,data1,markerl,markern,chr):
    chrdata = open(data1,'w')
    nsample = 0
    with open(cro, 'r') as cross:
        s1   = re.compile(r'^-s')
        s2   = re.compile(r'^\d+')
        sn   = re.compile(r'^-n\s+(\d+)')
        sp   = re.compile(r'^-p\s+(\d+)')
        tag1 = 0
        tag2 = 0
        num  = 0
        for line in cross:
            if sn.search(line):
                m = sn.search(line)
                nsample = m.groups(0)[0]
            if sp.search(line):
                m = sp.search(line)
                nloci = m.groups(0)[0]
                chrdata.write('Distinct = '+str(nsample)+'\n')
                chrdata.write('Gene = '+str(nsample)+'\n')
                chrdata.write('Loci = '+str(markern)+'\n')
                chrdata.write('I=1 %treat data as SNPs'+'\n')
                chrdata.write('K = 2 %a-allele model with Haplotype Alleles by 1,2'+'\n')
                chrdata.write('Positions of loci:'+'\n')
                chrdata.write(markerl+'\n')
                chrdata.write('Haplotypes'+'\n') 
            if s1.search(line):
                tag1 = 1
            if s2.search(line) and tag1 == 1:
                tag2 = 1
                num  = -1
            num += 1
            if tag1 == 1 and tag2 == 1:
                if num == chr:
                    #print num
                    #print line
                    line = re.sub(r'\s+',r'',line)
                    line = re.sub(r'0',r'1',line)
                    chrdata.write('         '+line+' 1\n')
    chrdata.write('#\n')
    chrdata.close()

def datafilehead(gmap,chr):
    marker = []
    with open(gmap, 'r') as gmap:
        s0   = re.compile(r'^-b MarkerNames')
        s1   = re.compile(r'^-e MarkerNames')
        s2   = re.compile(r'^'+str(chr)+r'\b')
        tag1 = 0
        tag2 = 0
        num  = 0
        for line in gmap:
            if s0.search(line):
                tag1 = 1
            if s1.search(line):
                tag1 = 0
            if s2.search(line) and tag1 == 1:
                line = line.rstrip()
                unit = line.split(' ')
                pos  = unit[2]
                pos  = pos[2:]
                marker.append(int(pos))
    markernum  = len(marker)
    markerline = ' '.join(map(str,marker))
    return (markerline,markernum)

def dict_add( dict0, key1, key2, value ):
   if dict0.has_key(key1):
       temp  = dict0.copy()
       dict1 = temp[key1]
       dict1.update({key2:value})
       temp[key1] = dict1
   else:
       temp  = dict0.copy()
       dict1 = {key2:value}
       temp.update({key1:{key2:value}})

   return temp


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
    for r in sorted (map (int, breakpoint.keys())):
        r = str(r)
        total = 0
        rilsbpdata=[]
        rilshead  =['RILs']
        rilsbpdata.append('GNN' + r)
        breakpointc = breakpoint[r]
        for c in sorted (breakpointc.keys()):
            chromosome = int (c)
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
    rilsbpfh.close()

     
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
   
    #Chr_breakpoint()
    RIL_breakpoint(binfile,rilbp)
    #GeneralSum()
     
'''     
    for chr in range(1,13):
        datafile1  = args.output + '.Chr' + str(chr) +'.Datafile'
        markerline, markernum = datafilehead(mapfile,chr)
        datafilehaplotype(crossfile,datafile1,markerline,markernum,chr) 
'''
if __name__ == '__main__':
    main()

