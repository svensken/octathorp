#! /bin/bash


for i in {1..30}
do
    cd $i
    /home/svensken/octathorp/score-thorp.py &
    echo $!
    cd -
done


#/home/svensken/octathorp/status-update.py &

