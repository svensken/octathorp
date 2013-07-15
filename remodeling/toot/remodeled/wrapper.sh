#!/bin/bash

for i in {6..10}
do
    sleep .5 # time to crtl-c out
    (
      name=3ANS_3GZJ_$i
      mkdir -p $name
      cd $name
      #echo $i | xargs -n 1 -P 10 ../remodel_simple.sh
      #( ../remodel.sh | tee $i.log ) &
      ../remodel.sh
      echo $! >> ../jobs.pid
      cd ..
    ) &
done 

