#!/usr/bin/env python

from rosetta import *
from structural_alignment import kabsch_alignment
import os, sys, time, datetime

# hell of an annoyance, this pymol messiness
sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/')
import pymol

pymol.finish_launching()
pymol.cmd.do('run ~/Desktop/PyRosetta/PyMOLPyRosettaServer.py')
time.sleep(2)

init()


def find_matching_domains( family_dir, orig_pdb_file, residue1, residue2, matchy, loggy ):
    
    ## calculate reference jump
    # load reference pose
    orig_pose = Pose()
    try:
        pose_from_pdb( orig_pose, orig_pdb_file )
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
    orig_pose_atom_tree = orig_pose.atom_tree()
    orig_rt = orig_pose_atom_tree.get_stub_transform( stub1, stub2 )


    # build pdb file list
    pdb_file_list = [item for item in os.listdir( family_dir ) if item[-4:] == '.pdb']

    n=0
    t0 = time.time()
    for cur_pdb_file_name in pdb_file_list:
        n=n+1
        t1 = time.time()
        try:
            cur_pose = Pose()
            pose_from_pdb( cur_pose, family_dir + cur_pdb_file_name )
            pose_load_result = '.'
        except PyRosettaException:
            pose_load_result = '*'
            #TODO remove 'continue', pass --ignore-unrecognized-rsd flag instead
            continue

        stub_list = []
        for middle_residue in range(2, cur_pose.total_residue()):
            # first stub's central residue is #2
            # last stub's central residue is #last-1
            s1a1 = AtomID(1, middle_residue - 1)
            s1a2 = AtomID(1, middle_residue)
            s1a3 = AtomID(1, middle_residue + 1)
            stub = StubID(s1a1, s1a2, s1a3)
            stub_list.append(stub)

        # populate pose's secstruct
        DsspMover().apply(cur_pose) # populate secstruct
        full_ss = cur_pose.secstruct()
        # infotize
        pdb_info = cur_pose.pdb_info()
        # atom-tree-ify
        cur_pose_atom_tree = cur_pose.atom_tree()
        
        # TODO un-globalize me
        global rmsd

        i = 1
        for first_stub in stub_list:
            for second_stub in stub_list[ i: ]: # start at subsequent stub (to first_stub)

                cur_rt = cur_pose_atom_tree.get_stub_transform(first_stub, second_stub)
                rmsd = distance( orig_rt, cur_rt )

                if rmsd < 5:
                    # TODO need better way to access residue numbers
                    cur_first_res  = first_stub.atom(2).rsd()
                    cur_second_res = second_stub.atom(2).rsd()
                    first_res_string  = pdb_info.pose2pdb( cur_first_res )
                    second_res_string = pdb_info.pose2pdb( cur_second_res )
                    cur_loop_length  = cur_second_res - cur_first_res
                    #cur_loop_length = cur_pose.residue(cur_first_res).polymeric_sequence_distance(cur_second_res)

                    # parameters to both check against and record to file
                    chunk_ss = full_ss[ cur_first_res : cur_second_res ]
                    this_ss_L    = chunk_ss.count('L') / float( len(chunk_ss) )
                    this_ss_E    = chunk_ss.count('E') / float( len(chunk_ss) )
                    this_ss_H    = chunk_ss.count('H') / float( len(chunk_ss) )

                    jump_length  = cur_rt.get_translation().length

                    
                    # tolerances for good matches
                    #TODO define at top of script, exec string (easily editable, and we can print tolerances into all_transforms file later)
                    q_tot_res       = 0  <  cur_pose.total_residue() < 500 
                    q_jump_length   = 0  <  jump_length              < 8
                    q_loop_length   = 10 <  cur_loop_length          < 80

                    if cur_pdb_file_name != orig_pdb_id and True: #q_jump_length and q_loop_length:

                        good_match = [    
                            cur_pdb_file_name,
                            cur_first_res, 
                            cur_second_res, 
                            cur_loop_length,  
                            jump_length,  
                            this_ss_L,    
                            this_ss_E,    
                            this_ss_H  ]
                        # OUPUT
                        for item in good_match:
                            matchy.write( str(item)+',' )
                        matchy.write( '\n' )
                        matchy.flush()
                        
                        print cur_pdb_file_name
                        print orig_pdb_file
                        # TODO make this part optional
                        construct_pose_from_matching_domains( orig_pose,            # pose
                                                              residue1,             # int
                                                              residue2,             # int
                                                              cur_pose,             # pose
                                                              cur_first_res,        # int
                                                              cur_second_res   )    # int

            i = i + 1


    loggy.write( '\nrun took ' + str(time.time()-t0) + ' seconds total\n')
    loggy.flush()




def construct_pose_from_matching_domains( host_pose,    # pose
                                          host_res1,    # int
                                          host_res2,    # int
                                          guest_pose,   # pose
                                          guest_res1,   # int
                                          guest_res2 ): # int

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
    
    # define residues for new pose
    our_Ns =  [host_pose.residue( r ).xyz('N') for r in range(1, host_res1+1)]
    our_Ns += [guest_pose.residue( r ).xyz('N') for r in range(guest_res1, guest_res2+1)]
    our_Ns += [host_pose.residue( r ).xyz('N') for r in range(host_res2, host_pose.total_residue())]

    # check clashiness
    n=1
    furthitudes = []
    for first_xyz in our_Ns:
        for second_xyz in our_Ns[ n: ]:
           furthitudes += [ first_xyz.distance( second_xyz ) ]
        n = n + 1
    loggy.write(str(furthitudes[:5]))
    loggy.flush()
    
    n_clashes = sum( f < 1 for f in furthitudes )
    # ^ looks like sum( [True, False, True, True] )
    # might be faster than saving a list and then evaluating

    print 'HOST: ', host_pose
    print 'GUEST: ', guest_pose
    print "####################"
    print 'rmsd ', rmsd
    print "CLASHES: ", n_clashes
    print 'guest_res1 ', guest_res1
    print 'guest_res2 ', guest_res2
    print "####################"

    # construct new pose
    if n_clashes < 10:
        pass

    pymover = PyMOL_Mover() 

    # i hate the trouble behind getting filename. stupid, stupid, stupid.
    # hope that pymol chooses this as the name
    host_pose_name = host_pose.pdb_info().name().split('/')[-1].split('.')[0]
    guest_pose_name = guest_pose.pdb_info().name().split('/')[-1].split('.')[0] # this one won't be a path, will it?

    pymover.apply(host_pose)
    time.sleep(.1)
    # prettify the visualization
    pymol.cmd.hide('everything')
    pymol.cmd.show('cartoon')
    time.sleep(.1)
    # easiest to see is pink cartoon for host pose and yellow lines for guest loop
    pymol.cmd.color('pink')

    pymover.apply(guest_pose)
    time.sleep(.1)
    pymol.cmd.hide('everything')
    pymol.cmd.show('cartoon')
    time.sleep(.1)
    pymol.cmd.color('yellow', guest_pose_name)

    # hide guest pose residues outside of the loop of interest
    pymol.cmd.hide( '( not resi ' + '+'.join( [str(t) for t in range(guest_res1-1, guest_res2+1)] ) + ' and '+guest_pose_name+' )' )
    # hide host pose residues outside of lame loop
    pymol.cmd.hide( '( '+host_pose_name+' and resi ' + '+'.join( [str(t) for t in range(host_res1, host_res2)] ) + ' )' )


    # generate new pose from aligned domains
    new_pose = Pose()


    raw_input('see pymol for match')
    
    # continue on to next match
    pymol.cmd.delete('all')



if __name__ == '__main__':

    with open( 'thorp-log', 'w') as loggy, open('matching_domains.'+str(time.strftime("%Y%m%d-%H%M%S")), 'w') as matchy:

        try:
            family_dir = sys.argv[-4]
            orig_pdb_file = sys.argv[-3]                # ./path/to/aaaa.pdb
            orig_pdb_id = orig_pdb_file.split('/')[-1]  # aaaa.pdb
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
            'invoked at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'+ \
            'Reference PDB: '+str(orig_pdb_file)+'\n'+ \
            'Residues '+str(residue1)+' and '+str(residue2)+'\n'+ \
            str(len(os.listdir(family_dir)))+' dir items (PDBs to posify and traverse for transforms)\n\n'
        loggy.write(output_string)
        loggy.flush()

        # TRANSFORMS OUTPUT
        # TODO output tolerances
        matchy_string = \
            'This list from '+family_dir+' directory, ### pdbs, *** tolerances, etc \n' + \
            'saved at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n'+ \
            str(sys.argv)+'\n\n\n'
        matchy.write(matchy_string)
        matchy.flush()
        
        find_matching_domains( family_dir, orig_pdb_file, residue1, residue2, matchy, loggy )

        #view_matching_domains()
        ## cycles through by default


