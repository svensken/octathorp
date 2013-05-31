#!/usr/bin/env python

print "importing rosetta"
from rosetta import *
import time, sys, os, csv

opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res', '-extra_res_cen', '/home/svensken/octathorp/molfiles/1WM1_PTB.cen.params', '-extra_res_fa', '/home/svensken/octathorp/molfiles/1WM1_PTB.fa.params' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )



this_directory = os.path.basename( os.getcwd() )

if not this_directory.isdigit():
    sys.exit(1)
    # dirtiest way to make sure we're in a dock output dir


# run from wrapper script

cr = csv.reader(open('energies.sc','rb'))

#              = [ (-63.3, 'manual_0.pdb), (...)]
all_cst_scores = [(float(row[0]),row[2]) for row in cr]
all_cst_scores.sort()
lowest_energy_filename = all_cst_scores[0][1]

rmsd_ref_pose = Pose()
pose_from_pdb( rmsd_ref_pose, os.path.basename(lowest_energy_filename) ) 

# why does cr empty itself upon first call?
cr = csv.reader(open('energies.sc','rb'))
for row in cr:

    cen_score_with_cst      = row[0]
    cen_score_without_cst   = row[1]
    filename                = row[2]

    this_pdb = os.path.join( this_directory, filename )

    pose = Pose()
    pose_from_pdb( pose, os.path.basename(filename) )

    # RMSD
    rmsd = CA_rmsd( rmsd_ref_pose, pose )

    with open('rmsds_energies.sc', 'a') as csv:
        linetowrite = str(rmsd)+','+cen_score_with_cst+','+cen_score_without_cst+','+filename
        csv.write( linetowrite + '\n' )



