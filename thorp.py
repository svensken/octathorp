#!/usr/bin/env python

from rosetta import *
import os, sys, time, numpy, datetime, operator
import pp

init()
#init(['app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ) ]) # pass flags here
# (pass flag to ignore unrecognized residues)


def find_similar_transforms( pdb_file, residue1, residue2 ):

    raw_input("please run 'tail -f loggy' to see all script output, then press enter")
    
    #TODO check chain of res
    #TODO disect workings of distance()
    #TODO pickle the overall jump_list to begin with
    #TODO detect domains (for res select without input)
    ## --> ExPASy's Prosite, Conserved Domain db, ESTHER db, ...

    ##################################################
    #- Get reference jump----------------------------#

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

    pose_atom_tree = pose.atom_tree()
    reference_jump = pose_atom_tree.get_stub_transform(stub1, stub2)

    #------------------------------------------------#
    ##################################################



    ##################################################
    #- Get all other jumps and compare --------------#

    # load all other PDB's
    # cleaned pdb's end with 'A'
    clean_pdb_file_list = [item for item in os.listdir('alpha-beta-hydrolases/') if item[-5:] == 'A.pdb']

    output_string = '\n'+str(len(clean_pdb_file_list))+' PDBs to posify and check jumps from:\n'
    f.write(output_string)
    f.flush()

    #clean_pose_list = []
    # to debug:
    global best_jumps
    best_jumps = []
    # load the poses one by one, check each against reference jump
    n=0
    t0 = time.time()
#TODO remove list limit
    for clean_pdb_file_name in clean_pdb_file_list[:3]: # first 100 files only
        n=n+1
        t1 = time.time()
        try:
            new_pose = Pose()
            pose_from_pdb( new_pose, 'alpha-beta-hydrolases/' + clean_pdb_file_name )
            pose_load_result = '.'
        except PyRosettaException:
            pose_load_result = '*'

        stub_list = []
        for middle_residue in range(2,new_pose.total_residue()):
            # first stub's central residue is #2
            # last stub's central residue is #last-1
            # Stub 1
            n_s1a1 = AtomID(1, middle_residue - 1)
            n_s1a2 = AtomID(1, middle_residue)
            n_s1a3 = AtomID(1, middle_residue + 1)
            n_stub = StubID(n_s1a1, n_s1a2, n_s1a3)
            stub_list.append(n_stub)

        new_pose_atom_tree = new_pose.atom_tree()
        i = 0 # or should be 1?
        global first_stub
        for first_stub in stub_list:
            for second_stub in stub_list[ i: ]: # start at subsequent stub (to first_stub)

                new_jump = new_pose_atom_tree.get_stub_transform(first_stub, second_stub)
                #new_jump = [ float(t) for t in str(new_jump).split() [1:] ] # omit 'RT'
                #big_diff = [ (a-b) / b for a,b in zip(new_jump, reference_jump) ]
                #avg_diff = numpy.average(big_diff)
                rms_deviation = distance(reference_jump, new_jump) 

                pdb_info = new_pose.pdb_info() # to store relevant info
                if rms_deviation < 20:
                    # structure of a good_jump:
                    # [ pdb, res1, res2, diff ]
                    ## what does rsd() return in terms of different chains?
                    ### should be pose's rsd, one sequence of numbers encapsulating all chains
                    # name could also be clean_pdb_file_name
                    good_jump = [                   \
                        pdb_info.name(),            \
                        pdb_info.pose2pdb( first_stub.atom(2).rsd() ),   \
                        pdb_info.pose2pdb( second_stub.atom(2).rsd() ),  \
                        rms_deviation ]
                    best_jumps.append(good_jump)
            i = i + 1

        if n%50==0 or n==len(clean_pdb_file_list):
            f.write('\n###'+str(n)+'###\n')
        else:
            f.write(str(round(time.time()-t1,3))+'s, ')
        f.flush()
        # this is where we need to pickle the Pose objects
        #clean_pose_list.append(pose)


    best_10_jumps_string = sorted( best_jumps, key=operator.itemgetter(3) )[ :10 ]

    output_string = \
        '\n\nAlles Klar\n\
         Time taken to load and grade all possible jumps: '+str(time.time()-t0)+' seconds\n'+ \
        'Best 10 jumps: \n\n'+ \
        '\n'.join( [str(L) for L in best_10_jumps_string] )
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

        find_similar_transforms( pdb_file, residue1, residue2 )
            
# close log file
f.close()

