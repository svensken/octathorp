#!/usr/bin/env python

import os, glob

combo_maker = {
    "1WM1" : [135, 98],
    "2WUF" : [133, 71],
    "3ANS" : [130, 105],
    "3B12" : [128, 87],
    "3GZJ" : [102, 74]  }

combo = "3ANS_3ANS"

# loopB
print combo
print "  bot:", combo[:4]
print "  cap:", combo[-4:]
botA = combo_maker[ combo[:4] ][0]
capB = combo_maker[ combo[:4] ][0] + 1 + combo_maker[ combo[-4:] ][1]
print "    botA:", botA
print "    capB:", capB

for length in ['5','6','7','8','9']:
    ampersand_command_B = "/home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
                    -database /home/svensken/Rosetta/main/database/ \
                    -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
                    -parser:script_vars anchorB="+ str(capB) +" \
                    -parser:protocol /home/svensken/octathorp/xmls/3ANS_by_length/loop_creation.B."+length+".xml \
                    -in:file::s /home/svensken/octathorp/xmls/3ANS_by_length/UM_4_D335H524D496Y465Y384_output_3ANS_1.pdb.massaged.pdbA.pdb"
    os.chdir('B.'+length)
    os.system( ampersand_command_B )#+ " &" )
    os.chdir('..')
    raw_input("length "+length+" complete")

raw_input("roundB done, check quantity (should be < 100*5)")

# loopA
print combo
print "  bot:", combo[:4]
print "  cap:", combo[-4:]
botA = combo_maker[ combo[:4] ][0]
capB = combo_maker[ combo[:4] ][0] + 1 + combo_maker[ combo[-4:] ][1]
print "    botA:", botA
print "    capB:", capB

count=0
for refined in glob.glob("*creation*refined*"):
    count += 1
    for length in ['8','9','10','11','12']:
        ampersand_command_A = "/home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
                        -database /home/svensken/Rosetta/main/database/ \
                        -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
                        -parser:script_vars anchorA="+ str(capA) +" \
                        -parser:protocol /home/svensken/octathorp/xmls/3ANS_by_length/loop_creation.A."+length+".xml \
                        -in:file::s " + refined
        if count >= 20:
            raw_input("A-count hit 20....")
            count=0
        os.chdir('A.'+length)
        os.system( ampersand_command_A )#+ " &" )
        oc.chdir('..')
        raw_input("length "+length+" complete for: "+refined)


# ls only fully completed structures
#pipe_to_xargs = "ls *creation*refined* | xargs -n1 -P20 "+command

