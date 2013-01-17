#!/usr/bin/env python

#from rosetta import Pose, pose_from_pdb
#from rosetta.core import kinematics
#from rosetta.core.kinematics import AtomTree
#from rosetta import init, Pose, pose_from_pdb, AtomID, StubID
from rosetta import *
import sys
init()


# input 'file1.pdb', res1, res2
try:
    pdb_file = [arg for arg in sys.argv if arg[-4:] == ".pdb"][0]
    print "accepted file: " + pdb_file
except:
    print "\nUsage: \n\
    $ thorp.py [file.pdb] [int residue 1] [int residue 2] \n\
    \n\
    Please put all relevant PDB files in the current directory.\n"
    #Please indicate path to PDB file relative to working directory
    sys.exit(0)


# load the PDB
pose = Pose()
try:
    pose_from_pdb( pose, 'sample_pdbs/' + pdb_file )
except:
    print "alles klar"

print pose

# recieve as input
residue1 = 22
residue2 = 33


#####
# select the atoms based on input
# Stub 1
s1a1 = AtomID(1, residue1 - 1)
s1a2 = AtomID(1, residue1)
s1a3 = AtomID(1, residue1 + 1)
stub1 = StubID(s1a1, s1a2, s1a3)
# Stub 2
s2a1 = AtomID(1, residue2 - 1)
s2a2 = AtomID(1, residue2)
s2a3 = AtomID(1, residue2 + 1)
stub2 = StubID(s2a1, s2a2, s2a3)

print stub1
print stub2


# how to find RT between stubs
# 1: RT core.kinematics.jump.jump(stub1,stub2), "ZERO rb_delta"
# 2: void core.kinematics.jump.from_stubs(stub1,stub2), "note: we don't reset rb_center!!!"
# 3: RT core.kinematics.RT.RT(stub1,stub2)
# 4: RT core.kinematics.AtomTree.get_stub_transform(stub1,stub2)
pose_atom_tree = pose.atom_tree()
jump = pose_atom_tree.get_stub_transform(stub1, stub2)

transl_list = [float(n) for n in str(jump.get_translation()).split()]
rotatn_list = [float(n) for n in str(jump.get_rotation()).split()]

# friggin all other RTs
# headless benchmark?


# check match against tolerance
transl_diff = [a/b for a,b in zip(transl_list1, transl_list2)]
rotatn_diff = [a/b for a,b in zip(rotatn_list1, rotatn_list2)]

comprehensive_diff = transl_diff + rotatn_diff
if [diff for diff in comprehensive_diff 

# would be faster with izip and immediate checks to break


########
# deem loop's secondary composition
# how to confine? just limit string to residue numbers?
ss = pose.secstruct()
percent_helical = 100. * ss.count("H") / len(ss)
#H=helical, E=sheet, L=loop
print str(percent_helical) + "% Helical"
