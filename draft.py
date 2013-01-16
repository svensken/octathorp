#!/usr/bin/env python

#from rosetta import Pose, pose_from_pdb
#from rosetta.core import kinematics
from rosetta import *
init()

# load the PDB
pose = Pose()
pose_from _pdb( pose, 'file.pdb' )

print pose.total_residue()
print pose.residue(33)
print pose.alpha(33)

# deem loop's secondary composition
# how to confine? just limit string to residue numbers?
ss = pose.secstruct()
percent_helical = 100. * ss.count("H") / len(ss)
#H=helical, E=sheet, L=loop

# fold tree
fold_tree = pose.fold_tree()
# why won't fold_tree.num_jump() return all possible?
# because a jump just connects chains


