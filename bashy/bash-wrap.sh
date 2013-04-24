#! /bin/bash


for i in {1..50}
do
    mkdir -p $i
    cd $i
    /home/svensken/octathorp/dock-thorp.py &
    echo $!
    cd -
done
