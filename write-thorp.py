#!/usr/bin/env python

from rosetta import *
import os, sys, time, numpy, datetime, operator
import pp, pickle

init()
#init(['app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ) ]) # pass flags here
# (pass flag to ignore unrecognized residues)


def get_all_transforms( ):

    raw_input("please run 'tail -f loggy' to see all script output, then press enter")
    
    #TODO detect domains (for res select without input)
    ## --> ExPASy's Prosite, Conserved Domain db, ESTHER db, ...


    ##################################################
    #- Get all other jumps and compare --------------#

    # load all other PDB's
    # cleaned pdb's end with 'A'
    pdb_file_list = [item for item in os.listdir('alpha-beta-hydrolases/') if item[-5:] == 'A.pdb']

    output_string = '\n'+str(len(pdb_file_list))+' PDBs to posify and check jumps from:\n'
    f.write(output_string)
    f.flush()

    best_jumps = []
    # load the poses one by one, check each against reference jump
    n=0
    t0 = time.time()
#TODO remove list limit
    for pdb_file_name in pdb_file_list[:3]: # first 100 files only
        n=n+1
        t1 = time.time()
        try:
            pose = Pose()
            pose_from_pdb( pose, 'alpha-beta-hydrolases/' + pdb_file_name )
            pose_load_result = '.'
        except PyRosettaException:
            pose_load_result = '*'

        stub_list = []
        for middle_residue in range(2, pose.total_residue()):
            # first stub's central residue is #2
            # last stub's central residue is #last-1
            # Stub 1
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
                this_ss_L    = chunk_ss.count('L') / len(chunk_ss)
                this_ss_E    = chunk_ss.count('E') / len(chunk_ss)
                this_ss_H    = chunk_ss.count('H') / len(chunk_ss)

                # define tolerances
                q_len  = pose.total_residue() < 500 # just to speed up first run
                # minimum E and H in ss
                
                if True and q_len:
                    # structure of a good_jump:
                    # [ jump, pdb-id, res1, res2, ss_L, ss_E, ss_H ]
                    good_jump = [           \
                        this_jump,          \
                        this_pdb_id,        \
                        this_1st_res,       \
                        this_2nd_res,       \
                        this_ss_L,          \
                        this_ss_E,          \
                        this_ss_H  ]
                    best_jumps.append(good_jump)
            i = i + 1

        if n%50==0 or n==len(pdb_file_list):
            f.write('\n###'+str(n)+'###\n')
        else:
            f.write(str(round(time.time()-t1,3))+'s, ')
        f.flush()
        # this is where we need to pickle the Pose objects
        #clean_pose_list.append(pose)


    best_10_jumps_string = sorted( best_jumps, key=operator.itemgetter(3) )[ :10 ]

    return best_10_jumps_string

    output_string = \
        '\n\nAlles Klar\n\
         Time taken to load and grade all possible jumps: '+str(time.time()-t0)+' seconds\n\n'+ \
        'Best 10 jumps: \n'+ \
        '' #'\n'.join( [str(L) for L in best_10_jumps_string] )
    f.write(output_string)
    f.flush()


        ## Alternatives;
        # check if 0 first
        # (max-min)/max
        # shift all +10
        # (need a better way)
        #print ''
        #print str(jamp[0]).split()[5]#10
        #print str(jamp[1]).split()[5]#10
        #print 'ref: ', reference_jump
        #print 'tmp: ', temp_jump
        #print 'diff list: ', big_diff
        #print 'avg diff:  ', avg_diff


    #------------------------------------------------#
    ##################################################


if __name__ == '__main__':

    with open( 'loggy', 'w') as f:
        
        try:
            # just some lazy defaults
            pdb_file = '1whtA.pdb' #[arg for arg in sys.argv if arg[-4:] == ".pdb"][0]
            residue1 = 22 #int(sys.argv[-2])
            residue2 = 33 #int(sys.argv[-1])
        except:
            print "\nUsage: \n\
            $ ./thorp.py [file.pdb] [int residue 1] [int residue 2] \n"
            sys.exit(0)
        
        # parallelize
        #job_server = pp.Server()
        #for cup in job_server.get_ncpus():

        best_jumps = get_all_transforms( )
            
# close log file
f.close()

