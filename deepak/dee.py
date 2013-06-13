#!/usr/bin/env python

from rosetta import *
from rosetta.protocols.rigid import *
import structural_alignment

opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-ignore_unrecognized_res', '-extra_res_fa', '/home/svensken/deepak/TST.params', '/home/svensken/deepak/AP1.params' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )


pose = Pose()
pose_from_pdb(pose, 'UM_1_S115H289D233_1EI9_clean_shdnew_2.pdb')

pep = Pose()
pose_from_pdb(pep, 'TST1_0001.pdb')

print pose.residue(280).atom('O1').xyz()
print pep.residue(1).atom('O9').xyz()



structural_alignment.kabsch_alignment( pep, pose, [ 15 , 20, 25 ], [ 10, 20, 30] )


print pose.residue(280).atom('O1').xyz()
print pep.residue(1).atom('O9').xyz()


dump_pdb pose
dump_pdb pep

system cat pepfile >> posefile

print 'good'
