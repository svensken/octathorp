#!/usr/bin/env python

import os

combo_maker = {
    "1WM1" : [135, 98],
    "2WUF" : [133, 71],
    "3ANS" : [130, 105],
    "3B12" : [128, 87],
    "3GZJ" : [102, 74]  }

os.chdir("combos/")

# loopB
for combo in os.listdir("."):
    if combo != "3ANS_3ANS": continue
    os.chdir( combo )

    print combo
    print "  bot:", combo[:4]
    print "  cap:", combo[-4:]
    botA = combo_maker[ combo[:4] ][0]
    capB = combo_maker[ combo[:4] ][0] + 1 + combo_maker[ combo[-4:] ][1]
    print "    botA:", botA
    print "    capB:", capB

    command = "/home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
                    -database /home/svensken/Rosetta/main/database/ \
                    -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
                    -parser:script_vars anchorB="+ str(capB) +" \
                    -parser:protocol /home/svensken/octathorp/xmls/loop_creation.B.xml \
                    -in:file::s "

    pipe_to_xargs = "ls | xargs -n1 -P20 "+command
    #os.system("ls")
    os.system( pipe_to_xargs )
    #raw_input("piped...")

    os.chdir('..')

raw_input("loopB done")

# loopA
for combo in os.listdir("."):
    if combo != "3ANS_3ANS": continue
    os.chdir( combo )

    print combo
    print "  bot:", combo[:4]
    print "  cap:", combo[-4:]
    botA = combo_maker[ combo[:4] ][0]
    capB = combo_maker[ combo[:4] ][0] + 1 + combo_maker[ combo[-4:] ][1]
    print "    botA:", botA
    print "    capB:", capB

    command = "/home/svensken/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease \
                    -database /home/svensken/Rosetta/main/database/ \
                    -lh:db_path /home/svensken/octathorp/xmls/2011vall/ \
                    -parser:script_vars anchorA="+ str(botA) +" \
                    -parser:protocol /home/svensken/octathorp/xmls/loop_creation.A.xml \
                    -in:file::s "

    # ls only fully completed structures
    pipe_to_xargs = "ls *creation*refined* | xargs -n1 -P20 "+command
    #os.system("ls")
    os.system( pipe_to_xargs )
    #raw_input("piped...")

    os.chdir('..')
