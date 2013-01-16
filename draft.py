#!/usr/bin/env python

#from rosetta import Pose, pose_from_pdb
#from rosetta.core import kinematics
import sys
from rosetta import init, Pose, pose_from_pdb, AtomID, StubID
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

print "Stub 1: " + stub1
print "Stub 2: " + stub2





########
# deem loop's secondary composition
# how to confine? just limit string to residue numbers?
ss = pose.secstruct()
percent_helical = 100. * ss.count("H") / len(ss)
#H=helical, E=sheet, L=loop


