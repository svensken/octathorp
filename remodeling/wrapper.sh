#!/bin/bash

for i in {1..10}
do
    mkdir -p $i
    cd $i
    ../remodel.sh &
    cd ..
done


