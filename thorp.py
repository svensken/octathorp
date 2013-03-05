#!/usr/bin/env python

from rosetta import *
import os, sys, time, numpy, datetime, operator

init()


with open( 'thorp-log', 'w') as loggy, open('all_transforms.'+str(time.strftime("%Y%m%d-%H%M%S")), 'w') as RTs:

    try:
        family_dir = sys.argv[-4]
        pdb_file = sys.argv[-3]
        residue1 = int(sys.argv[-2])
        residue2 = int(sys.argv[-1])
    except:
        print "input not recognized                         \n\
               Usage:                                       \n\
               $ ./read-thorp.py [directory/] [pdb_file] [res1] [res2]   \n"
        sys.exit(1)
   
    # build pdb list to check against
    pdb_file_list = [item for item in os.listdir( family_dir ) if item[-4:] == '.pdb']

    # DEBUG OUTPUT
    output_string = '\n\n#######################################\n'+ \
        'This file contains debugging output for thorp.py\n'+ \
        'invoked at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n' \
        'Reference PDB: '+str(pdb_file)+'\n'+ \
        'Residues '+str(residue1)+' and '+str(residue2)+'\n'+ \
        str(len(pdb_file_list))+' PDBs to posify and traverse for transforms\n\n'
    loggy.write(output_string)
    loggy.flush()

    # TRANSFORMS OUTPUT
    RTs_string = \
        'This list from *** directory, ### pdbs, *** tolerances, etc \n' + \
        'saved at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'
    RTs.write(RTs_string)
    RTs.flush()

    n=0
    t0 = time.time()
    for pdb_file_name in pdb_file_list:
        n=n+1
        t1 = time.time()
        try:
            pose = Pose()
            pose_from_pdb( pose, family_dir + pdb_file_name )
            pose_load_result = '.'
        except PyRosettaException:
            pose_load_result = '*'
            #TODO remove 'continue', pass --ignore-unrecognized-rsd flag instead
            continue

        stub_list = []
        for middle_residue in range(2, pose.total_residue()):
            # first stub's central residue is #2
            # last stub's central residue is #last-1
            s1a1 = AtomID(1, middle_residue - 1)
            s1a2 = AtomID(1, middle_residue)
            s1a3 = AtomID(1, middle_residue + 1)
            stub = StubID(s1a1, s1a2, s1a3)
            stub_list.append(stub)

        # populate pose's secstruct
        DsspMover().apply(pose) # populate secstruct
        full_ss = pose.secstruct()
        # infotize
        pdb_info = pose.pdb_info()
        # atom-tree-ify
        pose_atom_tree = pose.atom_tree()

        i = 1
        for first_stub in stub_list:
            for second_stub in stub_list[ i: ]: # start at subsequent stub (to first_stub)

                jump = pose_atom_tree.get_stub_transform(first_stub, second_stub)

                chunk_ss = full_ss[ first_stub.atom(2).rsd() : second_stub.atom(2).rsd() ]

                # define things to store
                this_jump    = str(jump)
                this_pdb_id  = pdb_file_name
                this_1st_res = pdb_info.pose2pdb( first_stub.atom(2).rsd() )
                this_2nd_res = pdb_info.pose2pdb( second_stub.atom(2).rsd() )
                this_length  = jump.get_translation().length
                this_ss_L    = chunk_ss.count('L') / float( len(chunk_ss) )
                this_ss_E    = chunk_ss.count('E') / float( len(chunk_ss) )
                this_ss_H    = chunk_ss.count('H') / float( len(chunk_ss) )

                #TODO define at top of script, exec string (so we can print tolerances into all_transforms file later)
                # define tolerances
                q_nrsd  = pose.total_residue() < 500 
                q_len   = this_length < 8
                q_loop  = int(this_2nd_res.split()[0]) - int(this_1st_res.split()[0]) > 10
                # minimum E and H in ss
                
                if True: # and q_len and q_loop:
                    # structure of a good_jump:
                    # [ jump, pdb-id, res1, res2, length, ss_L, ss_E, ss_H ]
                    good_jump = [           \
                        this_jump,          \
                        this_pdb_id,        \
                        this_1st_res,       \
                        this_2nd_res,       \
                        this_length,        \
                        this_ss_L,          \
                        this_ss_E,          \
                        this_ss_H  ]
                    # OUPUT
                    for item in good_jump:
                        RTs.write( str(item)+',' )
                    RTs.write( '\n' )
                    RTs.flush()
            i = i + 1

        if n%50==0 or n==len(pdb_file_list):
            loggy.write( '\n' + pdb_file_name )
            loggy.write('\n###'+str(n)+'###\n')
        #else:
        #    loggy.write(str(round(time.time()-t1,3))+'s, ')
        loggy.flush()

        
    
    
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


f = open('loggy-read-alt','w')

pymover = PyMOL_Mover() 

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
#with open('all_transforms', 'r') as giant_list:
with open('loggy-read', 'r') as giant_list:
    t0 = time.time()
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
                f.write( '###' + row_items[3] + ' ' + str(rmsd) + '\n' )
            if rmsd < 1:
                f.write(line + str(rmsd)+'\n\n')
                #construct_pose_from_matching_domains( pose,                # pose
                                                      residue1,             # int
                                                      residue2,             # int
                                                      row_list[1],          # '1aaaA.pdb'
                                                      int(row_items[2]),    # int
                                                      int(row_items[3])  )  # int"""

        except Exception,e:
            f.write(str(e)+'\n')
        if i % 1000000 == 0:
            f.write( str(i) + '\n')
        f.flush()
    
    # OUTPUT
    output_string = \
        '\n\nAlles Klar\n' + \
        'total run time: '+str(time.time()-t0)+' seconds\n\n'
    loggy.write(output_string)
    loggy.flush()

    f.write( 'read took ' + str(time.time()-t0) + ' seconds total')

f.close()
