#!/usr/bin/env python

from rosetta import *
import os, sys, time, datetime

init()


def find_matching_domains( family_dir, pdb_file, residue1, residue2, out_file=matchy, debug_file=loggy ):
    
    ## calculate reference jump
    # load reference pose
    pose = Pose()
    try:
        pose_from_pdb( pose, pdb_file )
        pose_load_result = 'no error'
    except PyRosettaException:
        pose_load_result = 'slight error'
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


    # build pdb file list
    pdb_file_list = [item for item in os.listdir( family_dir ) if item[-4:] == '.pdb']

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

                rt = pose_atom_tree.get_stub_transform(first_stub, second_stub)
                rmsd = distance( ref_rt, rt )

                if rmsd < 1:
                    # parameters to both check against and record to file
                    chunk_ss = full_ss[ first_stub.atom(2).rsd() : second_stub.atom(2).rsd() ]
                    this_ss_L    = chunk_ss.count('L') / float( len(chunk_ss) )
                    this_ss_E    = chunk_ss.count('E') / float( len(chunk_ss) )
                    this_ss_H    = chunk_ss.count('H') / float( len(chunk_ss) )

                    jump_length  = jump.get_translation().length

                    # TODO need better way to access residue numbers
                    this_1st_res = pdb_info.pose2pdb( first_stub.atom(2).rsd() )
                    this_2nd_res = pdb_info.pose2pdb( second_stub.atom(2).rsd() )
                    loop_length  = int(this_2nd_res.split()[0]) - int(this_1st_res.split()[0])

                    # tolerances for good matches
                    #TODO define at top of script, exec string (easily editable, and we can print tolerances into all_transforms file later)
                    q_tot_res       = pose.total_residue()      < 500 
                    q_jump_length   = jump_length               < 8
                    q_loop_length   = loop_length               < 80

                    if q_jump_length and q_loop_length:

                        good_match = [          \
                            pdb_file_name,      \
                            this_1st_res,       \
                            this_2nd_res,       \
                            jump_length,        \
                            loop_length,        \
                            this_ss_L,          \
                            this_ss_E,          \
                            this_ss_H  ]
                        # OUPUT
                        for item in good_match:
                            matchy.write( str(item)+',' )
                        matchy.write( '\n' )
                        matchy.flush()
                        
                        # TODO make this part optional
                        construct_pose_from_matching_domains( pose,                # pose
                                                              residue1,             # int
                                                              residue2,             # int
                                                              row_list[1],          # '1aaaA.pdb'
                                                              int(row_items[2]),    # int
                                                              int(row_items[3])  )  # int

            i = i + 1


    loggy.write( '\nrun took ' + str(time.time()-t0) + ' seconds total\n')
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
    kabsch_alignment( host_pose, guest_pose, [ host_res1 - 1 ,      \
                                               host_res1     ,      \
                                               host_res1 + 1 ,      \
                                               host_res2 - 1 ,      \
                                               host_res2     ,      \
                                               host_res2 + 1   ],   \
                                             [ guest_res1 - 1 ,     \
                                               guest_res1     ,     \
                                               guest_res1 + 1 ,     \
                                               guest_res2 - 1 ,     \
                                               guest_res2     ,     \
                                               guest_res2 + 1  ] )  \

    pymover = PyMOL_Mover() 
    pymover.apply(host_pose)
    pymover.apply(guest_pose)
    raw_input('see pymol for match')

    # generate new pose from aligned domains
    new_pose = Pose()

    f.write( str(host_pose)  )
    f.write( str(guest_pose) )
    f.flush()




f.close()

if __name__ == '__main__':

    with open( 'thorp-log', 'w') as loggy, open('matching_domains.'+str(time.strftime("%Y%m%d-%H%M%S")), 'w') as matchy:

        try:
            family_dir = sys.argv[-4]
            pdb_file = sys.argv[-3]
            residue1 = int(sys.argv[-2])
            residue2 = int(sys.argv[-1])
        except:
            print "input not recognized                                     \n\
                   Usage:                                                   \n\
                   $ ./thorp.py [directory/] [pdb_file] [res1] [res2]       \n\
                                                                            \n\
                   example: $ ./thorp.py super-family/ 3pf8A.pdb 138 182    \n"
            sys.exit(1)
       
        # DEBUG OUTPUT
        output_string = '\n\n#######################################\n'+ \
            'This file contains debugging output for thorp.py\n'+ \
            'invoked at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n' \
            'Reference PDB: '+str(pdb_file)+'\n'+ \
            'Residues '+str(residue1)+' and '+str(residue2)+'\n'+ \
            str(len(os.listdir(family_dir))+' PDBs to posify and traverse for transforms\n\n'
        loggy.write(output_string)
        loggy.flush()

        # TRANSFORMS OUTPUT
        # TODO output tolerances
        matchy_string = \
            'This list from '+family_dir+' directory, ### pdbs, *** tolerances, etc \n' + \
            'saved at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'
        matchy.write(matchy_string)
        matchy.flush()
        
        find_matching_domains( family_dir, pdb_file, residue1, residue2 ) # out_ and debug_files need not be passed

        #view_matching_domains()
        ## cycles through by default


