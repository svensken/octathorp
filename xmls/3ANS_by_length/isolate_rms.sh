#!/bin/bash


for pdb in UM*
do
    num=$(grep "rmsA " $pdb)
    echo $pdb","${num#* } >> rmseses
done

