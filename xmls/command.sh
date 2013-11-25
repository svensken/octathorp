#!/bin/bash 


# cutpoints
## ( name botA cap_len )
p1wm1=(1WM1 135 98)
p2wuf=(2WUF 133 71)
p3ans=(3ANS 130 105)
p3b12=(3B12 128 87)
p3gzj=(3GZJ 102 74)
combos=(p1wm1 p2wuf p3ans p3b12 p3gzj)


for combo_item in ${combos[@]}
do
    x="${combo_item}[@]"
    ints=( "${!x}" )

    echo name    ${ints[0]}
    echo botA    ${ints[1]}
    echo cap_len ${ints[2]}
    
    commandtorun="/home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
                        -database /home/svensken/Rosetta/main/database/ \
                        -lh:db_path 2011vall/ \
                        -parser:script_vars anchorA="+5+" anchorB="+7+" \ 
                        -parser:protocol loop_creation.xml \
                        -in:file::s "
    #ls combos/$combo | xargs $commandtorun
done

