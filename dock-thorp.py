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
hpose = Pose()
pose_from_pdb(hpose, '4fwb.pdb')
gpose = Pose()
pose_from_pdb(gpose, '2yas.pdb')
hres1 = 29
hres2 = 33
gres1 = 44
gres2 = 77

pi=hpose.pdb_info()

pymover = PyMOL_Mover()


# CHAINS
#TODO figure out chain manipulations to do this efficiently

#delete hres1-hres2
for r in reversed(range(hres1, hres2+1)): # WILL SEGFAULT AT r = 1
    hpose.delete_polymer_residue( r )

#append_residue_by_jump(hpose, gpose.res, hres1, start_new_chain=True)
hpose.append_residue_by_jump(gpose.residue(gres1), hres1, start_new_chain=True)

#append res of hpose reses
for r in reversed(range(gres1+1, gres2+1)):
    hpose.append_polymer_residue_after_seqpos(gpose.residue(r), hpose.total_residue(), 0)
    pi.chain( hpose.total_residue(), "B" )
    #pymover.apply(hpose)
    #time.sleep(.1)

hpose.dump_pdb('temp.pdb')
pose=Pose()
pose_from_pdb(pose, 'temp.pdb')

print pose.pdb_info()
print pose.fold_tree()
pymover.apply(pose)


# MOVER
randomize1 = RigidBodyRandomizeMover(pose)
#randomize2 = RigidBodyRandomizeMover(pose)

pert_mover = RigidBodyPerturbMover( 1, 8, 3 ) #(jump_num, rotation, translation)

spin = RigidBodySpinMover(1)
slide_into_contact = DockingSlideIntoContact(1
)

scorefxn_low = create_score_function('interchain_cen')
scorefxn_low.set_weight( core.scoring.atom_pair_constraint, 1 )

movemap = MoveMap()
movemap.set_jump(1, True)
minmover = MinMover()
minmover.movemap(movemap)
minmover.score_function(scorefxn_low)


perturb = SequenceMover()
perturb.add_mover(randomize1)
#perturb.add_mover(randomize2)
perturb.add_mover(pert_mover)
perturb.add_mover(spin)
perturb.add_mover(slide_into_contact)
perturb.add_mover(minmover)


####################
# in pose:
# hres1 is hres1
hterm1 = AtomID(1, hres1)
# hres2 is hres1+1
hterm2 = AtomID(1, hres1 + 1)
# gres1 is first B
gterm1 = AtomID(1, int(str(pose.fold_tree().jump_edge(1)).split()[2]) ) # please find a better way
# gres2 is last B
gterm2 = AtomID(1, pose.total_residue() )
####################

SOG = constraints.SOGFunc( 8, 1.0 )
swf = constraints.ScalarWeightedFunc( 1.0, SOG )

# gres1 and hres1
apc = constraints.AtomPairConstraint( hterm1, gterm1, swf )
pose.add_constraint( apc )

# gres2 and hres2
apc = constraints.AtomPairConstraint( hterm2, gterm2, swf )
pose.add_constraint( apc )


docking_low = DockingLowRes(scorefxn_low, 1)

#test_pose = Pose()
#test_pose.assign(pose)
AddPyMolObserver(pose, True)

# a. set necessary variables for this trajectory
# -reset the test pose to original (centroid) structure
# -change the pose name, for pretty output to PyMOL
pose.pdb_info().name('O_O')

# b. perturb the structure for this trajectory
print "now moving"
perturb.apply(pose)

print "now docking"
# c. perform docking
docking_low.apply(pose)
