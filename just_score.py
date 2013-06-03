#!/usr/bin/env python

import time, sys, os

print 'importing rosetta...'
from rosetta import *
from rosetta.protocols.rigid import *
# silence output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res', '-extra_res_cen', '/home/svensken/octathorp/molfiles/1WM1_PTB.cen.params', '-extra_res_fa', '/home/svensken/octathorp/molfiles/1WM1_PTB.fa.params' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )


pose = Pose()
pose_from_pdb( pose, sys.argv[1] )
# for benchmark, let's just define them manually
# POSE NUMBERING
hterm1 = 137 #140
hterm2 = 138 #251
gterm1 = 205
gterm2 = 310

hterm1 = AtomID( 1, hterm1 ) #hres1
hterm2 = AtomID( 1, hterm2 ) #hres2
gterm1 = AtomID( 1, gterm1 )#int(str(pose.fold_tree().jump_edge(1)).split()[2]) ) #chainB first res
gterm2 = AtomID( 1, gterm2 )#pose.total_residue()) #chainB last res

#to_centroid = SwitchResidueTypeSetMover('centroid')
#to_centroid.apply(pose)

#sf = create_score_function('score12')
sf = create_score_function_ws_patch('standard', 'score12')
sf.set_weight( atom_pair_constraint, 1 )

movemap = MoveMap()
movemap.set_bb(False)
#movemap.set_bb_true_range(pose.total_residue() / 4, pose.total_residue() * 3 / 4)
movemap.set_chi(False)
movemap.set_jump( 2, True )
#movemap.show(pose.total_residue())

minmover = MinMover()
minmover.movemap(movemap)
minmover.score_function(sf)

GF1 = constraints.GaussianFunc( 14.1, 2.0 )
GF2 = constraints.GaussianFunc( 11.5, 2.0 )
apc1 = constraints.AtomPairConstraint( hterm1, gterm1, GF1 )
apc2 = constraints.AtomPairConstraint( hterm2, gterm2, GF2 )
pose.add_constraint( apc1 )
pose.add_constraint( apc2 )


#pymover = PyMOL_Mover() 
#pymover.link.tcp_ip = "192.12.88.133" 
#pymover.link.tcp_port = 65000
#pymover.apply(pose)


sf.show(pose)

print 'before: ', sf(pose)

minmover.apply(pose)
pose.dump_pdb('minimized.pdb')

print 'after: ', sf(pose)

sf.show(pose)
