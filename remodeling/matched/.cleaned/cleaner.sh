#!/bin/bash

# run in remodeling/matched/.cleaned/

for f in $(ls ..)
do
    sed '/^$/d' ../$f > $f
    ../../toot/clean_pdb.py $f A
done


