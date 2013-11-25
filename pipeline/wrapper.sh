#!/bin/bash

combo_of_interest=$1
blueprint=$2

waitforme=()
for i in {1..10}
do
  remodel_call="\
/home/svensken/Rosetta/main/source/bin/remodel.default.linuxgccrelease \
-database /home/svensken/Rosetta/main/database/ \
-s /home/svensken/octathorp/pipeline/ready_to_remodel/"$combo_of_interest" \
-remodel:blueprint /home/svensken/octathorp/pipeline/ready_to_remodel/"$blueprint" \
-remodel:num_trajectory 1 \
-nstruct 2 \
-ex1 \
-ex2 \
-overwrite"
#-remodel:quick_and_dirty \

  thisdir="dir__"$combo_of_interest"__"$blueprint"__proc_"$i
  mkdir -p $thisdir
  cd $thisdir
  ( $remodel_call ) &
  waitforme+=($!)
  cd ..
done

wait ${waitforme[@]}

