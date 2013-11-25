#!/bin/bash


for length in 7 #{5..9}
do

    cd A.$length.noamp.rms
    /home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
        -database /home/svensken/Rosetta/main/database/ \
        -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
        -parser:protocol /home/svensken/octathorp/xmls/3ANS_by_length/loop_creation.A.xml \
        -parser:script_vars anchorA=130 lengthA=$length rmsEND=137 \
        -in:file::s /home/svensken/octathorp/xmls/3ANS_by_length/UM_4_D335H524D496Y465Y384_output_3ANS_1.pdb.massaged.pdbA.pdb #\
        #-in:file:native /home/svensken/octathorp/3ANS.pdbA.pdb
    cd ..

done

