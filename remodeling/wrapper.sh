#!/bin/bash

for i in {1..5}
do
    sleep .5 # time to crtl-c out
    (
      mkdir -p $i
      cd $i
      #echo $i | xargs -n 1 -P 10 ../remodel_simple.sh
      #( ../remodel.sh | tee $i.log ) &
      echo remodel #../remodel.sh
      echo $! >> ../jobs.pid
      cd ..
    ) &
done 

