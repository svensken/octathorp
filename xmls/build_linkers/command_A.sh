#!/bin/bash


for input in UM*creation*refined* 
do

    /home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
        -database /home/svensken/Rosetta/main/database/ \
        -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
        -parser:protocol /home/svensken/octathorp/xmls/build_linkers/loop_creation.A.xml \
        -parser:script_vars anchorA=130 lengthA=$1 tries=200 \
        -in:file:s $input
        #-in:file:s /home/svensken/octathorp/xmls/3ANS_by_length/UM_4_D335H524D496Y465Y384_output_3ANS_1.pdb.massaged.pdbA.pdb

done
