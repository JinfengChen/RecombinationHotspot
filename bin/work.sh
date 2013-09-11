python PreSequenceLDhot.py --input ../input/MPR.cross.uniq

echo "Summary map for RILs"
python MapSummary.py --input ../input/MPR.cross.uniq --bin ../input/MPR.geno.bin.uniq > HEG4svNB.log

echo "Hotspot using binomial test"
python SumHotspot.py --input ../input/MPR.geno.bin.uniq

