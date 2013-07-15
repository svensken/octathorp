#!/bin/bash

# run in remodeling/matched/.cleaned/

for f in $(ls ..)
do
    ../../toot/clean_pdb.py ../$f A
done


