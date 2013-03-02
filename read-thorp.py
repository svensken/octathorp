import sys, time
from rosetta import *
from structural_alignment import kabsch_alignment

init()

# TODO add support for multiple chains


def construct_pose_from_matching_domains( host_pose,    # pose
                                          host_res1,    # int
                                          host_res2,    # int
                                          guest_name,   # '1aaaA.pdb'
                                          guest_res1,   # int
                                          guest_res2 ): # int
    # build guest pose
    guest_pose = Pose()
    pose_from_pdb( guest_pose, 'alpha-beta-hydrolases/'+guest_name )

    # rotate guest pose to align with host's ref residues
    kabsch_alignment( host_pose, guest_pose, [ host_res1 - 1 ,
                                               host_res1     ,
                                               host_res1 + 1 ,
                                               host_res2 - 1 ,
                                               host_res2     ,
                                               host_res2 + 1   ],
                                             [ guest_res1 - 1 ,
                                               guest_res1     ,
                                               guest_res1 + 1 ,
                                               guest_res2 - 1 ,
                                               guest_res2     ,
                                               guest_res2 + 1  ] )
    pymover.apply(host_pose)
    pymover.apply(guest_pose)
    raw_input('see pymol for match')

    # generate new pose from aligned domains
    new_pose = Pose()

    f.write( str(host_pose)  )
    f.write( str(guest_pose) )
    f.flush()


f = open('loggy-read','w')

pymover = PyMOL_Mover() 

## GET MANUALLY-DEFINED REFERENCE JUMP
try:
    pass
    #pdb_file = [arg for arg in sys.argv if arg[-4:] == '.pdb'][0]
    #residue1 = int(sys.argv[-2]
    #residue2 = int(sys.argv[-1]
except:
    print "input not recognized                         \n\
           Usage:                                       \n\
           $ ./read-thorp.py [pdb_file] [res1] [res2]   \n"
    sys.exit(1)

# for now;
pdb_file = '3pf8A.pdb'
residue1 = 138
residue2 = 182

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
    'Residues '+str(residue1)+' and '+str(residue2)+'\n'+ \
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
# Transform
pose_atom_tree = pose.atom_tree()
ref_rt = pose_atom_tree.get_stub_transform( stub1, stub2 )


# reconstruct and check each RT object from giant list of jumps
with open('all_transforms', 'r') as giant_list:
    i=0
    for line in giant_list:
        i+=1
        try:
            row_items = line.split(',')
            rt_string = row_items[0].split()

            m = numeric.xyzMatrixdouble(0)
            m = m.rows( float(rt_string[1]),
                        float(rt_string[2]),
                        float(rt_string[3]),
                        float(rt_string[4]),
                        float(rt_string[5]),
                        float(rt_string[6]),
                        float(rt_string[7]),
                        float(rt_string[8]),
                        float(rt_string[9]) )
            v = numeric.xyzVector_double( float(rt_string[10]),
                                          float(rt_string[11]),
                                          float(rt_string[12]) )
            rt = RT()
            rt.set_rotation(m)
            rt.set_translation(v)
            
            rmsd = distance( ref_rt, rt )
            
            if row_items[1] == '3pf8A.pdb' and row_items[2] == '138 A':
                f.write( row_items[3] + ' ' + str(rmsd) + '\n' )
            if rmsd < 1:
                f.write(20*'BOOM')
                """construct_pose_from_matching_domains( pose,                # pose
                                                      residue1,             # int
                                                      residue2,             # int
                                                      row_list[1],          # '1aaaA.pdb'
                                                      int(row_items[2]),    # int
                                                      int(row_items[3])  )  # int"""

        except Exception,e:
            f.write(str(e)+'\n')
        if i % 1000 == 0:
            f.write( str(i) + '\n')
        f.flush()


f.close()
