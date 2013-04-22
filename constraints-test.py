#!/usr/bin/env python

print 'importing rosetta...'
from rosetta import *
from rosetta.protocols.rigid import *
import time


# silence rosetta output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )


# parameters
pose = Pose()
pose_from_pdb(pose, 'cst_test.pdb')
hterm1 = 136
hterm2 = 212
gterm1 = 217
gterm2 = 288

pi=pose.pdb_info()


#to_centroid = SwitchResidueTypeSetMover('centroid')
#to_centroid.apply(pose)


# Construct constraint set mover.
set_constraints = ConstraintSetMover()
set_constraints.constraint_file("temp.cst")
# Prepare scorefunction.
# one of: docking, docking_cen, docking_min, docking_occ_exact, empty, gauss, interchain_cen, 
sf = create_score_function("interchain_cen")
sf.set_weight(atom_pair_constraint, 10)
# Set constraints into pose.
set_constraints.apply(pose)
# Score the pose.
sf.show(pose)


#print sf









