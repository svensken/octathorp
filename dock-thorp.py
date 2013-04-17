#!/usr/bin/env python

print 'importing rosetta...'
from rosetta import *
import time


# silence rosetta's output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )


# parameters
hpose = Pose()
pose_from_pdb(hpose, '4fwb.pdb')
gpose = Pose()
pose_from_pdb(gpose, '2yas.pdb')
hres1 = 22
hres2 = 33
gres1 = 44
gres2 = 77

pi = hpose.pdb_info()
print pi


pymover = PyMOL_Mover()
#pymover.apply(hpose)
#raw_input('...')


# CHAINS
#delete hres1-hres2
#append_residue_by_jump(hpose, gpose.res, hres1, start_new_chain=True)
#append res of hpose reses

for r in reversed(range(hres1, hres2+1)): # WILL SEGFAULT AT r = 1
    hpose.delete_polymer_residue( r )

hpose.append_residue_by_jump(gpose.residue(gres1), hres1, start_new_chain=True)

for r in reversed(range(gres1+1, gres2+1)):
    hpose.append_polymer_residue_after_seqpos(gpose.residue(r), hpose.total_residue(), 0)
    pi.chain( hpose.total_residue(), "B" )
    pymover.apply(hpose)
    time.sleep(.1)

print pi


# EXPORT, IMPORT
hpose.dump_pdb('temp.pdb')
pose=Pose()
pose_from_pdb(pose, 'temp.pdb')

print pose.pdb_info()

pymover.apply(pose)



