python PreSequenceLDhot.py --input ../input/MPR.cross.uniq

echo "Summary map for RILs"
python MapSummary.py --input ../input/MPR.cross.uniq --bin ../input/MPR.geno.bin.uniq > HEG4svNB.log

echo "Hotspot using binomial test"
python SumHotspot.py --input ../input/MPR.cross.uniq --bin ../input/MPR.geno.bin.uniq

echo "LDhat"
python LDhat_SeqLDHot_Pipe.py --input ../input/BGI.SNP.Jap.matrix.1 --win 100000 > log 2> log2 &
python LDhat_SeqLDHot_Pipe.py --input ../input/BGI.SNP.Jap.matrix --win 100000 > log 2> log2 &
