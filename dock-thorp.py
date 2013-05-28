#!/usr/bin/env python

import time, os

print 'importing rosetta...'
from rosetta import *
from rosetta.protocols.rigid import *
# silence output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )

#print 'importing pymol...'
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/')
# silence output, no GUI
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
pymol.finish_launching()
#pymol.cmd.do('run ~/Desktop/PyRosetta/PyMOLPyRosettaServer.py')


def ce_align_poses():
    pymover = PyMOL_Mover()
    pymover.apply(hpose)
    pymover.apply(gpose)
    #time.sleep(1)
    pymol.cmd.load( hpose_filename )
    pymol.cmd.load( gpose_filename )
    hpose_name = '4fwb'
    gpose_name = '2yas'
    pymol.cmd.cealign(hpose_name, gpose_name)
    pymol.cmd.save('ce_host_temp.pdb', hpose_name)
    pymol.cmd.save('ce_guest_temp.pdb', gpose_name)
    #pose = Pose()
    #pose = Pose()
    pose_from_pdb(hpose, 'ce_host_temp.pdb')
    pose_from_pdb(gpose, 'ce_guest_temp.pdb')
    pymover.apply(hpose)
    pymover.apply(gpose)

    return hpose, gpose

def snip_poses():
    # delete between hres1 hres2
    # list should be [ hres2-1, hres2-2, ... , hres1+1 ]
    for r in reversed(range(hres1+1, hres2)): # WILL SEGFAULT AT r = 1
        hpose.delete_polymer_residue( r )

    # add guest domain to hpose as new chain
    #append_residue_by_jump(hpose, gpose.res, hres1, start_new_chain=True)
    hpose.append_residue_by_jump(gpose.residue(gres1), hres1, start_new_chain=True)
    pi=hpose.pdb_info()
    for r in reversed(range(gres1, gres2+1)):
        hpose.append_polymer_residue_after_seqpos(gpose.residue(r), hpose.total_residue(), 0)
        pi.chain( hpose.total_residue(), "B" )

    # reload to make new chain behave correctly
    hpose.dump_pdb('reload_temp.pdb')
    pose=Pose()
    pose_from_pdb(pose, 'reload_temp.pdb')

    return pose

def setup_movers_and_sf():
    pass

def setup_constraints():
    pass
    

##################################################

hpose_filename = '4fwb.pdb'
gpose_filename = '2yas.pdb'
hpose = Pose()
pose_from_pdb(hpose, hpose_filename)
gpose = Pose()
pose_from_pdb(gpose, gpose_filename)
hres1 = 136 #remains 136
hres2 = 212 #becomes 137
gres1 = 113 #becomes 219
gres2 = 185 #becomes 291


hpose, gpose = ce_align_poses()

pose = snip_poses()


sidechain_pose = Pose()
sidechain_pose.assign(pose)
to_centroid = SwitchResidueTypeSetMover('centroid')
to_centroid.apply(pose)


# SETUP MOVERS & SCOREFXN
## figure out how to wrap into a function
### (perturb would not return nicely... relies on the other variables in memory?)
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
#


# CONSTRAINTS
hterm1 = AtomID(1, hres1) #hres1
hterm2 = AtomID(1, hres1 + 1) #hres2
gterm1 = AtomID(1, int(str(pose.fold_tree().jump_edge(1)).split()[2]) ) #chainB first res
gterm2 = AtomID(1, pose.total_residue()) #chainB last res

GF = constraints.GaussianFunc( 6.0, 2.0 )
apc1 = constraints.AtomPairConstraint( hterm1, gterm1, GF )
apc2 = constraints.AtomPairConstraint( hterm2, gterm2, GF )
pose.add_constraint( apc1 )
pose.add_constraint( apc2 )

#cs = pose.constraint_set()
#print cs


starting_pose = Pose()
starting_pose.assign(pose)

print 'starting posed ok'

sf.show(pose)

print 'handshake ok'

# DOCKER & JOB-DISTRIBUTOR

docking_low = DockingLowRes(sf, 1)
docking_low.set_scorefxn( sf )

print 'scorefxn set ok'

AddPyMolObserver(pose, True)


#jd = PyJobDistributor('jd_output', 400, sf)
#jd.native_pose = pose

#while not jd.job_complete:
for a in range(200):
    print 'looped ok'
    # change pose name for PyMOL
    #pose.pdb_info().name('O_O')

    pose.assign(starting_pose) # has centroids and constraints

    print 'assigned ok'

    # perturb the structure
    perturb.apply(pose)

    print 'perturbed ok'

    # perform docking
    sf.show(pose)
    docking_low.apply(pose)
    sf.show(pose)

    recover_sidechains = ReturnSidechainMover(sidechain_pose)
    recover_sidechains.apply(pose)

    #jd.output_decoy(pose)

    # dump scored pdb (manually)
    filename = 'manual_' + str(a) + '.pdb'
    pose.dump_pdb( filename )#dump_scored_pdb( filename, sf )

    with open('status.update', 'a') as statusupdate:
        statusupdate.write(str(time.strftime("%Y/%m/%d-%H:%M:%S"))+'\n')
