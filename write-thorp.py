#!/usr/bin/env python

from rosetta import *
import os, sys, time, numpy, datetime, operator
import pp, pickle

init()
#init(['app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ) ]) # pass flags here
# (pass flag to ignore unrecognized residues)


def get_all_transforms( ):

    raw_input("please run 'tail -f loggy' to see clean script output, then press enter")
    
    # cleaned pdb's end with 'A'
    pdb_file_list = [item for item in os.listdir('alpha-beta-hydrolases/') if item[-5:] == 'A.pdb']

    # OUTPUT
    output_string = str(len(pdb_file_list))+' PDBs to posify and traverse for transforms\n\n'
    loggy.write(output_string)
    loggy.flush()

    n=0
    t0 = time.time()
    for pdb_file_name in pdb_file_list:
        n=n+1
        t1 = time.time()
        try:
            pose = Pose()
            pose_from_pdb( pose, 'alpha-beta-hydrolases/' + pdb_file_name )
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
                this_ss_L    = chunk_ss.count('L') #/ float( len(chunk_ss) )
                this_ss_E    = chunk_ss.count('E') #/ float( len(chunk_ss) )
                this_ss_H    = chunk_ss.count('H') #/ float( len(chunk_ss) )

                #TODO define at top of script, exec string (so we can print tolerances into all_transforms file later)
                # define tolerances
                q_nrsd  = pose.total_residue() < 500 
                q_len   = this_length < 8
                q_loop  = int(this_2nd_res.split()[0]) - int(this_1st_res.split()[0]) > 10
                # minimum E and H in ss
                
                if True and q_len and q_loop:
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
                    RTs.write( str([str(L) for L in good_jump]) )
                    RTs.flush()
            i = i + 1

        if n%50==0 or n==len(pdb_file_list):
            loggy.write('\n###'+str(n)+'###\n')
        else:
            loggy.write(str(round(time.time()-t1,3))+'s, ')
        loggy.flush()

        
    # OUTPUT
    output_string = \
        '\n\nAlles Klar\n' + \
        'total run time: '+str(time.time()-t0)+' seconds\n\n'
    loggy.write(output_string)
    loggy.flush()

    


if __name__ == '__main__':


    with open( 'loggy-write', 'w') as loggy, open('all_transforms', 'w') as RTs:
        #TODO backup RTs file if already populated
        
        # OUTPUT
        output_string = '\n\n#######################################\n'+ \
            'This file contains output from thorp.py\n'+ \
            'invoked at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'
        loggy.write(output_string)
        loggy.flush()

        # OUTPUT
        RTs_string = \
            'This list from *** directory, ### pdbs, *** tolerances, etc \n' + \
            'saved at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'
        RTs.write(RTs_string)
        RTs.flush()


        get_all_transforms()
            

