#!/usr/bin/env python

#from rosetta import Pose, pose_from_pdb
#from rosetta.core import kinematics
#from rosetta.core.kinematics import AtomTree
#from rosetta import init, Pose, pose_from_pdb, AtomID, StubID
from rosetta import *
# init, Pose, pose_from_pdb, kinematics, kinematics.AtomTree, AtomID, StubID
import sys


#sys.stdout = os.devnull
init() # assigning to variable won't silence it :(
#sys.stdout = sys.__stdout__



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
    silence = pose_from_pdb( pose, 'sample_pdbs/' + pdb_file )
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

ref_transl_list = [float(n) for n in str(jump.get_translation()).split()]
ref_rotatn_list = [float(n) for n in str(jump.get_rotation()).split()]


####################################################
#### Giant For Loop needed here ####################

# friggin all other RTs
# headless benchmark?

# test against self (try reverse)
new_pose = Pose()
try:
    pose_from_pdb(new_pose, 'sample_pdbs/' + pdb_file)
except:
    print 'good good'
total_residues = pose.total_residue()
new_pose_atom_tree = pose.atom_tree()

print 'clip here'

n = 0

for first_res_num in range(1,total_residues+1):
    for second_res_num in range(1,total_residues+1):
        # Stub 1
        n_s1a1 = AtomID(1, first_res_num - 1)
        n_s1a2 = AtomID(1, first_res_num)
        n_s1a3 = AtomID(1, first_res_num + 1)
        n_stub1 = StubID(n_s1a1, n_s1a2, n_s1a3)
        # Stub 2
        n_s2a1 = AtomID(1, second_res_num - 1)
        n_s2a2 = AtomID(1, second_res_num)
        n_s2a3 = AtomID(1, second_res_num + 1)
        n_stub2 = StubID(n_s2a1, n_s2a2, n_s2a3)
        
        n = n + 1
        print n
        
#new_jump = new_pose_atom_tree.get_stub_transform(n_stub1, n_stub2)


#temp_transl_list = [float(n) for n in str(jump.get_translation()).split()]
#temp_transl_list = [float(n) for n in str(jump.get_translation()).split()]

# check match against tolerance
#transl_diff = [a/b for a,b in zip(ref_transl_list, transl_list2)]
#rotatn_diff = [a/b for a,b in zip(ref_rotatn_list, rotatn_list2)]

#comprehensive_diff = transl_diff + rotatn_diff
#if [diff for diff in comprehensive_diff 
#for diff in comprehensive_diff:
#    if abs(diff - 1) >= .1 :
#        print "no match"
#        break
#    else:
#        print "good"
#        print stub1
#        print stub2

# would be faster with izip and immediate checks to break

########
# deem loop's secondary composition
# how to confine? just limit string to residue numbers?
#ss = pose.secstruct()
#percent_helical = 100. * ss.count("H") / len(ss)
#H=helical, E=sheet, L=loop
#print str(percent_helical) + "% Helical"

####################################################
#### For Loop done #################################

