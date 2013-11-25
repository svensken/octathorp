#!/bin/bash


for struct in "1WM1" "2WUF" "3ANS" "3B12" "3GZJ"
do
  for len in {3..8}
  do
    #/home/svensken/Rosetta/main/source/bin/score.default.linuxgccrelease -database /home/svensken/Rosetta/main/database -in:file:centroid -l lists/$struct.$len.list -in:file:native natives/$struct.pdb -score:weights cen_std -out:file:scorefile scores/$struct.$len.sc
    /home/svensken/Rosetta/main/source/bin/score.default.linuxgccrelease -database /home/svensken/Rosetta/main/database -in:file:l lists/$struct.$len.list -in:file:native natives/$struct.pdb -out:file:scorefile scores/$struct.$len.sc
  done
done
