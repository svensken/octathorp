#! /bin/bash


for i in {1..5}
do
    mkdir -p $i
    cd $!
    /home/svensken/octathorp/dock-thorp.py &
    cd -
done
