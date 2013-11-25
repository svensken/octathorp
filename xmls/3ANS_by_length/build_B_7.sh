#!/bin/bash


for length in {5..9}
do

cd B.$length.noamp
/home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
    -database /home/svensken/Rosetta/main/database/ \
    -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
    -parser:script_vars anchorB=236 \
    -parser:protocol /home/svensken/octathorp/xmls/3ANS_by_length/loop_creation.B.$length.xml \
    -in:file::s /home/svensken/octathorp/xmls/3ANS_by_length/UM_4_D335H524D496Y465Y384_output_3ANS_1.pdb.massaged.pdbA.pdb
cd ..

done
