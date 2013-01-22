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


# define all stubs THREE RESIDUES
#stub_list = []
#n = 1
#for middle_residue in range(2,new_pose.total_residue()):
#    # first stub's central residue is #2
#    # last stub's central residue is #last-1
#    print 'middle residue: ', middle_residue
#    n = n + 1
#    print 'iteration: ', n
#    # Stub 1
#    n_s1a1 = AtomID(1, middle_residue - 1)
#    n_s1a2 = AtomID(1, middle_residue)
#    n_s1a3 = AtomID(1, middle_residue + 1)
#    n_stub = StubID(n_s1a1, n_s1a2, n_s1a3)
#    print n_stub
#    stub_list.append(n_stub)

# define all stubs ONE RESIDUE
stub_list = []
for residue in range(1,new_pose.total_residue()+1):
    # first stub consists of atomno's 1, 2, and 3 of res#1
    # last stub consists of atomno's 1, 2, and 3 of res#last
    # Stub 1
    n_s1a1 = AtomID(1, residue - 1)
    n_s1a2 = AtomID(1, residue)
    n_s1a3 = AtomID(1, residue + 1)
    n_stub = StubID(n_s1a1, n_s1a2, n_s1a3)
    print n_stub
    stub_list.append(n_stub)


raw_input('hit enter to continue')
#print stub_list[57].atom(2) # middle


new_pose_atom_tree = new_pose.atom_tree()
jump_list = []
for first_stub in stub_list:
    for second_stub in stub_list:
        new_jump = new_pose_atom_tree.get_stub_transform(first_stub, second_stub)
        print new_jump
        jump_list.append(new_jump)


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

#find_stub_transform (?)
# last but not least;
# http://graylab.jhu.edu/pyrosetta/downloads/scripts/demos/D120_Ligand_interface.py
