#!/usr/bin/env python

print 'importing rosetta...'
from rosetta import *
from rosetta.protocols.rigid import *
import time, os

if False: #try:
    viz_ = True
    print 'importing pymol...'
    sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/')
    import __main__
    __main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
    import pymol
    pymol.finish_launching()
    pymol.cmd.do('run ~/Desktop/PyRosetta/PyMOLPyRosettaServer.py')

    pymol.cmd.do('print "heyllo"')

    #pymol.cmd.ray()
    #pymol.cmd.png
else: #except:
    viz_ = False
    print 'hmm, no pymol.'


# silence rosetta output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )

# path to script. seems like a dumb way to locate the pdbs.
hier = str(os.path.dirname(os.path.realpath(__file__))) # unicode weirdness?

# parameters
hpose = Pose()
pose_from_pdb(hpose, hier+'/4fwb.pdb')
gpose = Pose()
pose_from_pdb(gpose, hier+'/2yas.pdb')
hres1 = 136 #remains 136
hres2 = 212 #becomes 137
gres1 = 113 #becomes 219
gres2 = 185 #becomes 291


if viz_:
    pymover = PyMOL_Mover()
    hpose.pdb_info().name('O_O')
    pymover.apply(hpose)
    pymover.apply(gpose)
    raw_input('ready to delete hreses')



### POSES

pi=hpose.pdb_info()

# delete between hres1 hres2
# list should be [ hres2-1, hres2-2, ... , hres1+1 ]
for r in reversed(range(hres1+1, hres2)): # WILL SEGFAULT AT r = 1
    hpose.delete_polymer_residue( r )
    if viz_:
        hpose.pdb_info().name('O_O')
        pymover.apply(hpose)

if viz_:
    raw_input('ready to add greses as chain B')

#append_residue_by_jump(hpose, gpose.res, hres1, start_new_chain=True)
hpose.append_residue_by_jump(gpose.residue(gres1), hres1, start_new_chain=True)

#append res of hpose reses
for r in reversed(range(gres1, gres2+1)):
    hpose.append_polymer_residue_after_seqpos(gpose.residue(r), hpose.total_residue(), 0)
    pi.chain( hpose.total_residue(), "B" )
    if viz_:
        hpose.pdb_info().name('O_O')
        pymover.apply(hpose)


hpose.dump_pdb('temp.pdb')
pose=Pose()
pose_from_pdb(pose, 'temp.pdb')

# (before centroid)
sidechain_pose = Pose()
sidechain_pose.assign(pose)

to_centroid = SwitchResidueTypeSetMover('centroid')
to_centroid.apply(pose)

if viz_:
    pose.pdb_info().name('O_O')
    pymover.apply(pose)




### MOVERS & SCOREFXN

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


### CONSTRAINTS

hterm1 = AtomID(1, hres1) #hres1
hterm2 = AtomID(1, hres1 + 1) #hres2
gterm1 = AtomID(1, int(str(pose.fold_tree().jump_edge(1)).split()[2]) ) #chainB first res
gterm2 = AtomID(1, pose.total_residue()) #chainB last res

#better way to get hterms
# assume chB goes res r to total_residue()
#r = 1
#while pi.chain(r) != 'B':
#r += 1
#hterm1_res = pose.residue( hres1 )                  #136
#hterm2_res = pose.residue( hres1+1 )                #137
#
#gterm1_res = pose.residue( r+1 )                    #219
#gterm2_res = pose.residue( pose.total_residue() )   #291


#swf = constraints.ScalarWeightedFunc( 10, SOG ) #what dis be?
GF = constraints.GaussianFunc( 6.0, 2.0 )
apc1 = constraints.AtomPairConstraint( hterm1, gterm1, GF )
apc2 = constraints.AtomPairConstraint( hterm2, gterm2, GF )
pose.add_constraint( apc1 )
pose.add_constraint( apc2 )

#cs = pose.constraint_set()
#print cs

starting_pose = Pose()
starting_pose.assign(pose)
print_if_desired = sf.show(pose)


# DOCKER & JOB-DISTRIBUTOR

docking_low = DockingLowRes(sf, 1)
docking_low.set_scorefxn( sf )

AddPyMolObserver(pose, True)


jd = PyJobDistributor('jd_output', 400, sf)
jd.native_pose = pose

#while not jd.job_complete:
for a in range(400):
    # change pose name for PyMOL
    pose.pdb_info().name('O_O')

    pose.assign(starting_pose) # has centroids and constraints

    # perturb the structure
    print "now moving"
    perturb.apply(pose)

    # perform docking
    print "now docking"
    sf.show(pose)
    docking_low.apply(pose)

    sf.show(pose)

    jd.output_decoy(pose)

    recover_sidechains = ReturnSidechainMover(sidechain_pose)
    recover_sidechains.apply(pose)

    # dump scored pdb (manually)
    filename = 'manual_' + str(a) + '.pdb'
    pose.dump_scored_pdb( filename, sf )

    with open('status.update', 'a') as statusupdate:
        statusupdate.write('klar med rund '+str(a)+', kl.'+str(time.strftime("%Y%m%d-%H%M%S"))+'\n')
