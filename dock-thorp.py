#!/usr/bin/env python

print 'importing rosetta...'
from rosetta import *
from rosetta.protocols.rigid import *
import time, os


# silence rosetta output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )

hier = os.path.dirname(os.path.realpath(__file__))

# parameters
hpose = Pose()
pose_from_pdb(hpose, hier+'/4fwb.pdb')
gpose = Pose()
pose_from_pdb(gpose, hier+'/2yas.pdb')
hres1 = 136 #remains 136
hres2 = 212 #becomes 137
gres1 = 113 #becomes 219
gres2 = 185 #becomes 291

pi=hpose.pdb_info()


# CHAINS
#TODO find a better way to manipulate the chains

# delete between hres1 hres2
# list should be [ hres2-1, hres2-2, ... , hres1+1 ]
for r in reversed(range(hres1+1, hres2)): # WILL SEGFAULT AT r = 1
    hpose.delete_polymer_residue( r )

#append_residue_by_jump(hpose, gpose.res, hres1, start_new_chain=True)
hpose.append_residue_by_jump(gpose.residue(gres1), hres1, start_new_chain=True)

#append res of hpose reses
for r in reversed(range(gres1, gres2+1)):
    hpose.append_polymer_residue_after_seqpos(gpose.residue(r), hpose.total_residue(), 0)
    pi.chain( hpose.total_residue(), "B" )
    #pymover.apply(hpose)
    #time.sleep(.1)

hpose.dump_pdb(hier+'/temp.pdb')
pose=Pose()
pose_from_pdb(pose, hier+'/temp.pdb')
to_centroid = SwitchResidueTypeSetMover('centroid')
to_centroid.apply(pose)

print pose.pdb_info()
print pose.fold_tree()


# MOVER
randomize1 = RigidBodyRandomizeMover(pose)
randomize2 = RigidBodyRandomizeMover(pose)

pert_mover = RigidBodyPerturbMover( 1, 30, 20 ) #(jump_num, rotation, translation)

spin = RigidBodySpinMover(1)
slide_into_contact = DockingSlideIntoContact(1)

sf = create_score_function('interchain_cen')
sf.set_weight( atom_pair_constraint, 10 )


movemap = MoveMap()
movemap.set_jump(1, True)
minmover = MinMover()
minmover.movemap(movemap)
minmover.score_function(sf)


perturb = SequenceMover()
perturb.add_mover(randomize1)
perturb.add_mover(randomize2)
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
gterm2 = AtomID(1, pose.total_residue())
####################

# cst file
#with open('temp.cst', 'w') as cst_file:
#    cst_file.write( 'AtomPair CA '+str(hres1+1)+' CA '+str(pose.total_residue())+' GAUSSIANFUNC 8.0 2.0')
    #AtomPair: Atom1_Name Atom1_ResNum Atom2_Name Atom2_ResNum Func_Type Func_Def
    #Distance between Atom1 and Atom2
    #GAUSSIANFUNC: mean sd

# pyrosetta tutorial method
#sc = ConstraintSetMover()
#sc.constraint_file('temp.cst')
#sc.apply(pose)

# whats this?
#swf = constraints.ScalarWeightedFunc( 10, SOG )

GF = constraints.GaussianFunc( 6.0, 2.0 )
apc1 = constraints.AtomPairConstraint( hterm1, gterm1, GF )
apc2 = constraints.AtomPairConstraint( hterm2, gterm2, GF )
pose.add_constraint( apc1 )
pose.add_constraint( apc2 )


# show scores!
print sf.show(pose)


docking_low = DockingLowRes(sf, 1)
docking_low.set_scorefxn( sf )

AddPyMolObserver(pose, True)


starting_pose = Pose()
starting_pose.assign(pose)

jd = PyJobDistributor('output', 1, sf)
jd.native_pose = pose 

#jd.job_complete = False
print jd.job_complete

while not jd.job_complete:
#for a in range(20):
    # change pose name for PyMOL
    pose.pdb_info().name('O_O')

    pose.assign(starting_pose)

    # perturb the structure
    print "now moving"
    perturb.apply(pose)

    # perform docking
    print "now docking"
    docking_low.apply(pose)

    jd.output_decoy(pose)

    # dump scored pdb (manually)
#TODO make job distributor dump pdbs and scorefiles
    filename = os.getcwd() + '/scored_' + str(a) + '.pdb'
    pose.dump_scored_pdb( filename, sf )

    print "past last"
