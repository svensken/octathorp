#!/usr/bin/env python

print 'importing rosetta...'
from rosetta import *


# silence rosetta's output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )


# parameters
host_pose = Pose()
pose_from_pdb(host_pose, '')
guest_pose = Pose()
pose_from_pdb(guest_pose, '')
hres1 = 
hres2 = 
gres1 = 
gres2 = 


# CHAINS
#delete hres1-hres2
#append_residue_by_jump(hpose, gpose.res, hres1, start_new_chain=True)
#append res of hpose reses

