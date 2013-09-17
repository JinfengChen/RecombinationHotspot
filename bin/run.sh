#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=5gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

#/rhome/cjinfeng/software/tools/sequenceLDhot/sequenceLDhot Infile1 SeqLDhot.Chr12.Datafile SeqLDhot.Chr12.output
python LDhat_SeqLDHot_Pipe.py --input ../input/BGI.SNP.Jap.matrix --output BGI_Jap_multithread --win 1000000

echo "Done"
