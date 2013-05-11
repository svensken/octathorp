#!/usr/bin/env python

import os, json
from datetime import datetime


#TODO get nstruct and ncpu *programmatically*
#  ncpu can be read from number of dirs in bashy/
#  nstruct can be read from source code (??)
#  these are lame. i want integrated



# assume this script is run as soon as dock run starts
initial_timestamp = datetime.now()


def take_snapshot():
    
    json_data = []

    for root, dirs, files in os.walk("bashy-old"):
        
        dirname = os.path.basename(root) # '12'
        # skip unintended dirs
        if not dirname.isdigit():
            continue
        print dirname

        # number of decoys done
        ndone = len( [f for f in files if f.startswith('output') and f.endswith('.pdb')] )
        print ndone

        # fully worthless decoys
        losers = 0
        if False:#for filename in directory:
            with open( filename, 'r' ) as f:
                for line in f:
                    if '-nan' in line:
                        losers += 1
        # should also check for 
        
        # time predictions
        with open( os.path.join(root, 'status.update'), 'r') as timely:
            times = timely.readlines()
            delta_list = []
            for l in range(times):
                try:
                    t1 = datetime.strptime( times[l], "%Y/%m/%d-%H:%M:%S" )
                    t2 = datetime.strptime( times[l+1], "%Y/%m/%d-%H:%M:%S" )
                    delta = t2 - t1
                    delta_list.append( delta.seconds )
                except:
                    pass
            avg_delta = sum(delta_list) / len(delta_list)
            # could potentially be more statusupdates than pdbs (append mode in dock-thorp.py)
            npdbs_1hr = 3600 / avg_delta
            npdbs_10hr = 36000 / avg_delta
        print avg_delta


        # add a bar
        json_data.append( {
                            "title"     : str(dirname),
                            "subtitle"  : "decoys",
                            "ranges"    : "50",
                            "measures"  : "70",
                            "markers"   : "90"
                          } )


        with open("/home/svensken/octathorp/d3tst/otherbullets.json", 'w') as json_file:
            print 'a'#json.dump( json_data, json_file )


        print 'klar :D'





if True: # watchdog control goes here
    take_snapshot()
