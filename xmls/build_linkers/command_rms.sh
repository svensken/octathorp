#!/bin/bash

for pdb in $(cat /home/svensken/octathorp/xmls/3ANS_by_length/A_B_built.list); do
    /home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
        -database /home/svensken/Rosetta/main/database/ \
        -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
        -parser:protocol /home/svensken/octathorp/xmls/3ANS_by_length/rms.xml \
        -in:file:s $pdb \
        -in:file:native /home/svensken/octathorp/3ANS.pdbA.pdb
    done

