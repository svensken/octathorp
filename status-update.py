#!/usr/bin/env python

import json
from datetime import datetime


#TODO get nstruct and ncpu *programmatically*
#  ncpu can be read from number of dirs in bashy/
#  nstruct can be read from source code (??)
#  these are lame. i want integrated



# assume this script is run as soon as dock run starts
initial_timestamp = datetime.now()


def take_snapshot:
    
    json_data = []

    for directory in bashy:
        
        # finished decoys
        ndone = len( [a for a in os.listdir(directory) if a.startswith('manual') and a.endswith('.pdb')] )

        # fully worthless decoys
        losers = 0
        for filename in directory:
            with open( filename, 'r' ) as f:
                for line in f:
                    if '-nan' in line:
                        losers += 1
        # should also check for 
        
        # time predictions
        with open('status.update', 'r') as timely:
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



        # add a bar
        json_data.append( {
                            "title"     : str(dirname),
                            "subtitle"  : "decoys",
                            "ranges"    : hmm,
                            "measures"  : hmm,
                            "markers"   : hmm
                          } )




    with open("/home/svensken/octathorp/d3tst/bullets.json", 'w') as json_file:
        json.dump( json_data, json_file )


    print 'klar :D'
