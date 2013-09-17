#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import re
import os
import argparse
from numpy import *
sys.path.append("/rhome/cjinfeng/software/ProgramPython/module/lib")
from multi import worker

def usage():
    test="name"
    message='''
python LDhat_SeqLDHot_Pipe.py --input ../input/BGI.SNP.Jap.matrix --win 100000
Use genotype data from one population to estimate recombination rate in 100 kb windows. And Find recombination hotspot using sequenceLDhot in these 100 kbwindows.

Genotype data are input as matrix:
chromosome01 411 A A A A A A A A A A A A A A A A A A A G A A A
chromosome01 465 G G G G G G G G G G G G G G G G G G G A G G G
chromosome01 482 T T T T T T T T T T T T T T T T T T T C T T T

    '''
    print message

def matrix2LDhat(matrix, win,output):
    samplen = 0
    siten   = 0
    phase   = 1
    tag     = 0
    index   = 0
    chro    = ''
    locs    = []
    data    = defaultdict(list)
    s       = re.compile(r'(\d+)')
    subfile = []
    if not os.path.isdir(output):
        os.system('mkdir ' + output)
    with open (matrix, 'r') as matrixfh:
        for line in matrixfh:
            line = line.rstrip()
            unit = re.split(r'\s+',line)
            samplen = len(unit) - 2
            index   = int(unit[1])/int(win)
            m       = s.search(unit[0])
            chrn    = m.groups(0)[0]
            chro    = 'Chr' + str(chrn)
            #print unit[1], index
            snps    = unit[2:]
            for i in range(len(snps)):
                if index == tag:
                    data['Genotype'+str(i)].append(snps[i])
                else:
                    '''site file'''
                    outfile = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win)) + '.sites'
                    print outfile
                    prefix  = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win))
                    subfile.append(prefix) 
                    with open (outfile, 'w') as filefh:
                        print >> filefh, ' ' + str(samplen) + ' ' + str(siten) + ' ' + str(phase)
                        for g in sorted(data.keys()):
                            seq = ''.join(data[g])
                            print >> filefh, '>' + g
                            print >> filefh, seq
                    '''loc file'''
                    locfile = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win)) + '.locs'
                    print locfile
                    with open (locfile, 'w') as filefh:
                        print >> filefh, str(siten) + ' ' + str(win) + ' L'
                        locsline = ' '.join(locs)
                        print >> filefh, locsline
                    '''datafile for sequenceLDhot'''
                    datafile = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win)) + '.Datafile'
                    print datafile
                    with open (datafile, 'w') as filefh:
                        print >> filefh, 'Distinct = ' + str(samplen)
                        print >> filefh, 'Gene = ' + str(samplen)
                        print >> filefh, 'Loci = ' + str(siten)
                        print >> filefh, 'I = 1 %treat data as SNPs'
                        print >> filefh, 'K = -4 %4-allele model with Haplotype Alleles specified by A,C,G,T'
                        print >> filefh, 'Position of loci:'
                        newlocs = map (lambda x: str(int(float(x)*1000.00)), locs)
                        locsline = ' '.join(newlocs)
                        print >> filefh, locsline
                        print >> filefh, 'Haplotypes'
                        for g in sorted(data.keys()):
                            seq = ''.join(data[g])
                            print >> filefh, '         ' + seq + ' 1'
                        print >> filefh, '#'
                    '''clear parameter'''
                    tag   = index
                    siten = 0
                    locs  = []
                    data.clear()
                    siten += 1
                    data['Genotype'+str(i)].append(snps[i])
            locs.append(str(float(unit[1])/1000.00))
            siten = len(locs)
        '''output the last windows'''
        '''site file'''
        outfile = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win)) + '.sites'
        print outfile
        prefix  = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win))
        subfile.append(prefix)
        with open (outfile, 'w') as filefh:
            print >> filefh, ' ' + str(samplen) + ' ' + str(siten) + ' ' + str(phase)
            for g in data.keys():
                seq = ''.join(data[g])
                print >> filefh, '>' + g
                print >> filefh, seq
        '''site file'''
        locfile = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win)) + '.locs'
        print locfile 
        with open (locfile, 'w') as filefh:
            print >> filefh, str(siten) + ' ' + str(win) + ' L'
            locsline = ' '.join(locs)
            print >> filefh, locsline
        '''datafile for sequenceLDhot'''
        datafile = output + '/' + output + '.' + chro + '.' + str(tag * int(win)) + '_' + str((tag + 1) * int(win)) + '.Datafile'
        print datafile
        with open (datafile, 'w') as filefh:
            print >> filefh, 'Distinct = ' + str(samplen)
            print >> filefh, 'Gene = ' + str(samplen)
            print >> filefh, 'Loci = ' + str(siten)
            print >> filefh, 'I = 1 %treat data as SNPs'
            print >> filefh, 'K = -4 %4-allele model with Haplotype Alleles specified by A,C,G,T'
            print >> filefh, 'Position of loci:'
            newlocs = map (lambda x: str(int(float(x)*1000.00)), locs)
            locsline = ' '.join(newlocs)
            print >> filefh, locsline
            print >> filefh, 'Haplotypes'
            for g in sorted(data.keys()):
                seq = ''.join(data[g])
                print >> filefh, '         ' + seq + ' 1'
            print >> filefh, '#'
    return subfile

def runLDhat(subfile):
    #LDhat = '/Users/jinfengchen/biocluster/RecombinationHotspot/tools/LDhat_v2.2'
    LDhat = '/rhome/cjinfeng/software/tools/LDhat_v2.2'
    cmds  = []
    for prefix in subfile: 
        cmdline = LDhat + '/rhomap -seq ' + prefix + '.sites -loc ' + prefix + '.locs -lk lk23.txt -its 1000000 -samp 2000 -burn 100000 -exact -prefix ' + prefix + '.'
        print cmdline
        cmds.append(cmdline)
        #os.system(cmdline)                
    worker(cmds, 20)

def getregionrho(sumfile):
    rhos = []
    with open (sumfile, 'r') as filefh:
        header = filefh.readline()
        for line in filefh:
            line = line.rstrip()
            unit = re.split(r'\t',line)
            rhos.append(float(unit[2]))
    rho = "%.2f" %mean(rhos)
    return rho

def Infile(rho, infile1):
    text = '''Number of runs = 5000
MIN Number of iterations per hotspot = 100
driving values (for rho) = 2
background rho = ''' + str(rho) + '''
theta (per site) = 0.001
abs grid for hotspot likelihood
0.5 40
rel grid for hotspot likelihood
10 100
sub-region (number of SNPS; length (bps); frequency (bps))
7 2000 1000
#
'''
    with open (infile1, 'w') as filefh:
        print >> filefh, text

def sequenceLDhot(subfile):
    #LDhot = '/Users/jinfengchen/biocluster/RecombinationHotspot/tools/sequenceLDhot'
    LDhot = '/rhome/cjinfeng/software/tools/sequenceLDhot'
    cmds  = []
    for prefix in subfile:
        sumfile = prefix + '.summary.txt'
        infile1 = prefix + '.Infile1'
        datafile = prefix + '.Datafile'
        rho = getregionrho(sumfile)
        Infile(rho, infile1)
        cmdline = LDhot + '/sequenceLDhot ' + infile1 + ' ' + datafile
        print cmdline
        cmds.append(cmdline)
        #os.system(cmdline)
    worker(cmds, 20)

def summaryLDhot(prefix):
    Rcmd ='''
#summary sequenceLDhot
source("/rhome/cjinfeng/software/tools/sequenceLDhot/HotspotSummary.R")
HotspotSummary("''' + prefix + '/' + prefix  + '.*.sum", 10, 4, T, T, "' + prefix + '.LDhot.summary", "' + prefix + '''")

''' 
    Rscript = prefix+'.LDhot.R'
    with open (Rscript, 'w') as filefh:
        print >> filefh, Rcmd
    cmdline = 'cat ' + Rscript + ' | R --slave'
    os.system(cmdline)
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-w', '--win')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    if args.win is None:
        args.win = 100000
    if args.output is None:
        args.output = 'BGI_Jap'
 
    '''Convert the genotype from matrix to site and loc for LDhat and to Datafile for sequenceLDhot'''
    subfile = matrix2LDhat(args.input,args.win,args.output)
    runLDhat(subfile)
    sequenceLDhot(subfile)
    summaryLDhot(args.output) 
    

if __name__ == '__main__':
    main()

