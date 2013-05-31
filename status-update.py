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

    #TODO find way to extract nstruct
    nstruct = 300

    # for combo runs, jump into newest directory
    all_subdirs = ['production/'+d for d in os.listdir('production/') if os.path.isdir('production/'+d)]
    latest_subdir = max(all_subdirs, key=os.path.getmtime)
    print latest_subdir

    total_pdbs = 0
    total_scored = 0
    for root, dirs, files in os.walk( latest_subdir ):
        
        dirname = os.path.basename(root) # '12'
        
        # skip unintended dirs
        if not dirname.isdigit():
            continue
        #print 'd',dirname

        # number of decoys done
        ndone = len( [f for f in files if f.startswith('manual') and f.endswith('.pdb')] )
        #print 'n',ndone

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
            for l in range(len(times)):
                try:
                    t1 = datetime.strptime( times[l][-15:-1], "%Y%m%d-%H%M%S" )
                    t2 = datetime.strptime( times[l+1][-15:-1], "%Y%m%d-%H%M%S" )
                    delta = t2 - t1
                    delta_list.append( delta.seconds )
                except:
                    pass
            try:
                avg_delta = sum(delta_list) / len(delta_list)
                # could potentially be more statusupdates than pdbs (append mode in dock-thorp.py)
                n1hr = 3600 / avg_delta
                n10hr = 36000 / avg_delta
            except:
                avg_delta = 0
                n1hr = 0
                n10hr = 0
        #print 'a',avg_delta

        nneg = 20 # how to get number of negative energy pdbs?

        # number scored
        with open( os.path.join(root, 'rmsd_vs_energy.csv'), 'r') as scorelist:
            scores = scorelist.readlines()
            wc = len(scores) - 1 # number of lines

        # add a bar
        json_data.append( {
                            "title"     : latest_subdir + dirname,
                            "subtitle"  : "decoys",
                            "ranges"    : [wc,wc], #[n1hr,n10hr],
                            "measures"  : [ndone,ndone],
                            "markers"   : [nstruct]
                          } )

        # i'll plug this in somewhere
        total_pdbs += ndone
        total_scored += wc


    with open("/home/svensken/octathorp/d3tst/otherbullets.json", 'w') as json_file:
        json.dump( json_data, json_file )


    print 'klar :D'
    print 'total pdbs:', total_pdbs
    print 'total scored:', total_scored
    time_elapsed = datetime.now() - initial_timestamp
    print 'time elapsed:', time_elapsed

    # things to put in status page header thing
    # (bulletchart legend)
    # full_total
    # initial_timestamp and time_elapsed




if True: # watchdog control goes here
    take_snapshot()
