#!/opt/Python/2.7.3/bin/python
import sys
import re
import os
import argparse

def usage():
    message='''
python PreSequenceLDhot.py --input ../input/MPR.cross.uniq  --output pipe.conf
The script read *.map and *.cro file as input, parse the input file and output Datafile for SequenceLDhot.
--input: input is a prefix of .cro and .map. ../input/MPR.cross.uniq will stand for ../input/MPR.cross.uniq.cro and ../input/MPR.cross.uniq.map
--output: output is then prefix of output file. Rice will make the output file to Rice.Infile1 and Rice.Datafile. Optional, Defaule is SeqLDhot.

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
    if (args.output is None):
        args.output = 'SeqLDhot'       

    print 'INPUT file is' + args.input
    print 'OUTPUT file is' + args.output
    crossfile = args.input + '.cro'
    mapfile   = args.input + '.map'
    infile1   = args.output + '.Infile1'
    #datafile1  = args.output + '.Datafile'
    for chr in range(1,13):
        datafile1  = args.output + '.Chr' + str(chr) +'.Datafile'
        markerline, markernum = datafilehead(mapfile,chr)
        datafilehaplotype(crossfile,datafile1,markerline,markernum,chr) 

if __name__ == '__main__':
    main()

