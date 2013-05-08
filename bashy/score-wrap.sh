#! /bin/bash


for i in {1..10}
do
    cd $i
    /home/svensken/octathorp/score-thorp.py &
    echo $!
    cd -
done


/home/svensken/octathorp/status-update.py &

