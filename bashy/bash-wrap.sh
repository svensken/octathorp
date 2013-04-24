#! /bin/bash


for i in {1..5}
do
    mkdir -p $i
    cd $i
    /home/svensken/octathorp/dock-thorp.py &
    echo $!
    cd -
done
