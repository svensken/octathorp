#!/bin/bash


for struct in "1WM1" "2WUF" "3ANS" "3B12" "3GZJ"
do
for len in {3..8}
do
for i in {1..10}
do
  for a in "0001" "0002"
  do
    echo remodeled/dir__$struct.trim.pdbA.pdb__$struct.$len.blueprint__proc_$i/$struct.trim.pdbA_$a.pdb >> lists/$struct.$len.list
  done
done
done
done
