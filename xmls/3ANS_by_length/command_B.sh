#!/bin/bash


for length in 10 #{8..12}
do

    cd B.$length.A_built.200_B_each
    /home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
        -database /home/svensken/Rosetta/main/database/ \
        -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
        -parser:protocol /home/svensken/octathorp/xmls/3ANS_by_length/loop_creation.B.xml \
        -parser:script_vars anchorB=243 lengthB=$length \
        -in:file:l /home/svensken/octathorp/xmls/3ANS_by_length/A_built.list
    cd ..

done

