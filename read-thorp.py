from rosetta import *
import sys
init()


## GET MANUALLY-DEFINED REFERENCE JUMP
try:
    pdb_file = [arg for arg in sys.argv if arg[-4:] == '.pdb'][0]
    residue1 = int(sys.argv[-2]
    residue2 = int(sys.argv[-1]
except:
    print "input not recognized                         \n\
           Usage:                                       \n\
           $ ./read-thorp.py [pdb_file] [res1] [res2]   \n"
    sys.exit(1)

## load reference pose
pose = Pose()
try:
    pose_from_pdb( pose, 'alpha-beta-hydrolases/' + pdb_file )
    pose_load_result = 'no error'
except PyRosettaException:
    pose_load_result = 'slight error'

output_string = '\n\n##################################################\n'+ \
    'This file contains output from the thorp.py script\n'+ \
    'invoked at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'+ \
    'Reference PDB: '+str(pdb_file)+'\n'+ \
    'Reference jump residues: '+str(residue1)+' and '+str(residue2)+'\n'+ \
    str(pose_load_result)+' loading the reference pose (no need to worry)\n'+ \
    'total reference residues: '+str(pose.total_residue())+'\n\n'
f.write(output_string)
f.flush()

## calculate reference jump
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

#------------------------------------------------#
##################################################



# reconstruct an RT object from a stored string
m = numeric.xyzMatrixdouble(1)
m = m.rows(11,12,13,14,15,16,17,18,19)
v = numeric.xyzVector_double(1,2,3)

rt = RT()
rt.set_rotation(m)
rt.set_translation(v)

#print rt
#print rt.get_translation()
#print rt.get_rotation()



# compare
rmsd = distance( reference_jump, this_jump )
