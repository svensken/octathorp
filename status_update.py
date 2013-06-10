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
    nstruct = 1120 + 790

    # for combo runs, jump into newest directory
    #all_subdirs = ['testing/'+d for d in os.listdir('testing/') if os.path.isdir('testing/'+d)]
    #latest_subdir = max(all_subdirs, key=os.path.getmtime)
    #print latest_subdir

    total_pdbs = 0
    total_scored = 0
    if True: #for root, dirs, files in os.walk( latest_subdir ):
        
        #dirname = os.path.basename(root) # '12'
        
        # skip unintended dirs
        #if not dirname.isdigit():
        #    continue
        #print 'd',dirname

        # number of decoys done
        #ndone = len( [f for f in files if f.startswith('manual') and f.endswith('.pdb')] )
        with open('matching/all_matched/status.update','r') as a:
            ndone = int(a.read())
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
        if False: #with open( os.path.join(root, 'status.update'), 'r') as timely:
            times = timely.readlines()
            delta_list = []
        #try:
            for l in range(len(times)-1):
                t1 = datetime.strptime( times[l][:19], "%Y/%m/%d-%H:%M:%S" )
                t2 = datetime.strptime( times[l+1][:19], "%Y/%m/%d-%H:%M:%S" )
                delta = t2 - t1
                delta_list.append( delta.seconds )
            avg_delta = sum(delta_list) / len(delta_list)
            n5min = 300 / avg_delta
            n10min = 600 / avg_delta
            n1hr = 3600 / avg_delta
            n10hr = 36000 / avg_delta
            #print 'a '+str(avg_delta)+'s'

        #except:
            #avg_delta, n1hr, n10hr = 0,0,0
            # could potentially be more statusupdates than pdbs (append mode in dock-thorp.py)

        nneg = 20 # how to get number of negative energy pdbs?

        # number scored
        #try:
        #    with open( os.path.join(root, 'rmsd_vs_energy.csv'), 'r') as scorelist:
        #        scores = scorelist.readlines()
        #        wc = len(scores) - 1 # number of lines
        #except:
        #    wc = 0

        # add a bar
        json_data.append( {
                            "title"     : 'decoys', #latest_subdir+'/'+dirname,
                            "subtitle"  : "decoys",
                            "ranges"    : [nneg],
                            "measures"  : [ndone,ndone],
                            "markers"   : [nstruct]
                          } )

        # i'll plug this in somewhere
        total_pdbs += ndone
        #total_scored += wc


    with open("/home/svensken/octathorp/d3tst/otherbullets.json", 'w') as json_file:
        json.dump( json_data, json_file )


    print 'klar :D'
    print 'total pdbs:', total_pdbs
    #print 'total scored:', total_scored
    time_elapsed = datetime.now() - initial_timestamp
    print 'time elapsed:', time_elapsed

    # things to put in status page header thing
    # (bulletchart legend)
    # full_total
    # initial_timestamp and time_elapsed




if True: # watchdog control goes here
    take_snapshot()
