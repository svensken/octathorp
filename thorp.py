#!/usr/bin/env python

# TODO
#   - define selections of residues in pymol, makes subsequent manipulation much easier
#   - define tolerances at top of script, exec string (easily editable, and we can print tolerances into all_transforms file later)

print 'importing rosetta...'
from rosetta import *
from structural_alignment import kabsch_alignment
import os, sys, time, datetime
import argparse



# silence rosetta's output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts ) 
core.init( args )



def gimme_stubs( r1, r2, interval_len=1 ):
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
    
    given_stub_list = [ stub1, stub2 ]
    return given_stub_list

#def find_diverge( 


def find_matching_domains( args ):
    
    ## calculate reference jump
    # load reference pose
    orig_pose = Pose()
    try:
        pose_from_pdb( orig_pose, orig_pdb_file )
    except PyRosettaException:
        pass

    # get a flag call in here
    # for blindly checking all jumps in area
    given_stub_list = gimme_stubs( residue1, residue2 )
    stub1 = given_stub_list[0]
    stub2 = given_stub_list[1]
    
    # Transform
    orig_pose_atom_tree = orig_pose.atom_tree()
    orig_rt = orig_pose_atom_tree.get_stub_transform( stub1, stub2 )


    # build pdb file list
    if args.directory:
        pdb_file_list = [item for item in os.listdir( args.directory ) if item[-4:] == '.pdb']
    elif args.guest_pdb:
        pdb_file_list = [args.guest_pdb]
    else:
        pdb_file_list = []

    
    #######################################
    n=0
    t0 = time.time()
    for cur_pdb_file_name in pdb_file_list:
        n=n+1
        t1 = time.time()
        try:
            cur_pose = Pose()
            if args.directory:
                pose_from_pdb( cur_pose, args.directory + cur_pdb_file_name )
            elif args.guest_pdb:
                cur_pose = Pose()
                pose_from_pdb( cur_pose, cur_pdb_file_name )
        except PyRosettaException:
            # what pyrosettaexceptions exist other than unrecognized residues?
            continue
        

        # define residues to check
        if args.guest_pdb and args.res1 and args.res2:
            residue_check_list = [ args.res1, args.res2 ]
        else:
            residue_check_list = range(2, cur_pose.total_residue())
            # first stub's central residue is #2
            # last stub's central residue is #last-1

        stub_list = []
        for middle_residue in residue_check_list:
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
        
        # TODO unglobalize
        global rmsd

        i = 1
        for first_stub in stub_list:
            for second_stub in stub_list[ i: ]: # start at subsequent stub (to first_stub)

                cur_rt = cur_pose_atom_tree.get_stub_transform(first_stub, second_stub)
                rmsd = distance( orig_rt, cur_rt )

                # TOLERANCE (put at top)
                q_rmsd = rmsd < 5

                if rmsd < args.rmsd:

                    # TODO need better way to access residue numbers
                    cur_first_res  = first_stub.atom(2).rsd()
                    cur_second_res = second_stub.atom(2).rsd()

                    cur_loop_length  = cur_second_res - cur_first_res
                    #cur_loop_length = cur_pose.residue(cur_first_res).polymeric_sequence_distance(cur_second_res)
                    first_res_string  = pdb_info.pose2pdb( cur_first_res )
                    second_res_string = pdb_info.pose2pdb( cur_second_res )

                    # TOLERANCE (put at top)
                    chunk_ss = full_ss[ cur_first_res : cur_second_res ]
                    this_ss_L    = chunk_ss.count('L') / float( len(chunk_ss) )
                    #this_ss_E    = chunk_ss.count('E') / float( len(chunk_ss) )
                    #this_ss_H    = chunk_ss.count('H') / float( len(chunk_ss) )

                    #jump_length  = cur_rt.get_translation().length

                    # tolerances for good matches
                    q_loop_length = cur_loop_length < args.length
                    q_ss_L = this_ss_L * 10 < args.loop_percentage
                    #q_jump_length   = 0  <  jump_length              < 8

                    if cur_pdb_file_name != orig_pdb_id and q_loop_length and q_ss_L:

                        construct_pose_from_matching_domains( orig_pose,            # pose
                                                              residue1,             # int
                                                              residue2,             # int
                                                              cur_pose,             # pose
                                                              cur_first_res,        # int
                                                              cur_second_res   )    # int
                        #good_match = [    
                        #   cur_pdb_file_name,
                        #   cur_first_res, 
                        #   cur_second_res, 
                        #   cur_loop_length,  
                        #   jump_length,  
                        #   this_ss_L,    
                        #   this_ss_E,    
                        #   this_ss_H  ]
                        # OUPUT
                        #for item in good_match:
                        #    matchy.write( str(item)+',' )
                        #matchy.write( '\n' )
                        #matchy.flush()

            i = i + 1


    #loggy.write( '\nrun took ' + str(time.time()-t0) + ' seconds total\n')
    #loggy.flush()




def construct_pose_from_matching_domains( old_host_pose,    # pose
                                          host_res1,    # int
                                          host_res2,    # int
                                          old_guest_pose,   # pose
                                          guest_res1,   # int
                                          guest_res2 ): # int


    # copy poses
    # TODO is this step actually necessary? any significant speed toll?
    host_pose = Pose()
    host_pose.assign( old_host_pose )
    guest_pose = Pose()
    guest_pose.assign( old_guest_pose )


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
    host_rsds =  [r for r in range(1, host_res1+1)] + [r for r in range(host_res2, host_pose.total_residue())]
    guest_rsds = [r for r in range(guest_res1, guest_res2+1)]


    # check clashiness
    n=1
    close_ones = 0
    for gn in guest_rsds:
        for hn in host_rsds:
            dist = guest_pose.residue( gn ).xyz('CA').distance( host_pose.residue( hn ).xyz('CA') )
            if dist < 3.0:
                close_ones += 1
    
    if close_ones > .01*args.clash_percentage*len(guest_rsds):
        return 1
 

    # too much trouble just to get pdb_id.pdb, annoying.
    # hope that pymol chooses this as the name
    host_pose_name = host_pose.pdb_info().name().split('/')[-1].split('.')[0]
    guest_pose_name = guest_pose.pdb_info().name().split('/')[-1].split('.')[0]

    print 'HOST: ', host_pose_name
    print 'GUEST: ', guest_pose_name
    print "####################"
    print 'rmsd ', rmsd
    print "CLASHES: ", close_ones
    print 'guest_res1 ', guest_res1
    print 'guest_res2 ', guest_res2
    print "####################"


    if args.visualize:

        pymover = PyMOL_Mover()

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
        pymol.cmd.hide( '( not resi ' + str(guest_res1)+'-'+str(guest_res2) + ' and '+guest_pose_name+' )' )
        # hide host pose residues outside of lame loop
        pymol.cmd.hide( '( '+host_pose_name+' and resi ' + str(host_res1)+'-'+str(host_res2) + ' )' )

        # show old loop slightly transparent and gray?

        raw_input('hit enter for CE')
        # CEALIGN
        pymol.cmd.cealign(host_pose_name, guest_pose_name)
        # analyze
        pymol.cmd.save('ce_host.pdb', host_pose_name)
        #time.sleep?
        pymol.cmd.save('ce_guest.pdb', guest_pose_name)
        ce_host_pose = Pose()
        ce_guest_pose = Pose()
        pose_from_pdb(ce_host_pose, 'ce_host.pdb')
        pose_from_pdb(ce_guest_pose, 'ce_guest.pdb')

        print ce_host_pose
        print ce_guest_pose

        raw_input('hit enter to continue to next match')



    if args.outfiles:
        # already working with duplicate poses, just modify in place
        # to remove residues, delete backwards; rosetta updates sequence index to keep continuity from 1
        for r in reversed(range(host_res1, host_res2+1)): # WILL SEGFAULT AT r = 1
            host_pose.delete_polymer_residue( r )
        for r in reversed(guest_rsds):
            host_pose.append_polymer_residue_after_seqpos( guest_pose.residue( r ), host_res1, 0 )

        lego_pose = Pose()
        lego_pose.assign( host_pose )

        # TODO shorten pose name after debugged and stuff
        ref_pose = host_pose_name+'-'+str(host_res1)+'-'+str(host_res2)
        origin_pose = host_pose_name+'-'+str(host_res1)+'-'+str(host_res2)
        match_pose = guest_pose_name+'-'+str(guest_res1)+'-'+str(guest_res2)
        #match_params = '%.4f-%i' % (rmsd, close_ones)

        try:
            os.mkdir( 'pose-dumps/' )
        except OSError:
            pass # dir exists
        dump_name = 'pose-dumps/'+origin_pose+'_'+match_pose+'.pdb'
        lego_pose.dump_pdb( dump_name )

        print 'successfully dumped to file ' + dump_name







if __name__ == '__main__':

    
    parser = argparse.ArgumentParser()

    # not optional
    parser.add_argument( "ref_pdb", help="path to host pdb file" )
    parser.add_argument( "ref_res_1", type=int, help="residue 1 for reference pdb" )
    parser.add_argument( "ref_res_2", type=int, help="residue 2 for reference pdb" )

    m_group = parser.add_mutually_exclusive_group()
    # directory of guests
    m_group.add_argument( "-d", "--directory", help="directory of pdbs to check against" )
    # specific guest
    m_group.add_argument( "-g", "--guest_pdb", help="path to a single pdb to check against" )
    # specific reses in guest
    parser.add_argument( "--res1", type=int, help="use only these residues in guest_pdb" )
    parser.add_argument( "--res2", type=int, help="use only these residues in guest_pdb" )

    # output file
    parser.add_argument( "-o", "--outfiles", help="whether or not pdbs are created (in pose-dumps/ dir) from the matches", action="store_true" )
    # visualize
    parser.add_argument( "-V", "--visualize", help="use import-capable pymol (not MacPyMOL) to visualize matches as they are found", action="store_true" )
    # tolerances
    parser.add_argument( "--rmsd",              default=3,  help="maximum root-mean-squared deviation between transforms (default is 3)" )
    parser.add_argument( "--length",            default=100,help="maximum number of residues in matching loop (default is 100)" )
    parser.add_argument( "--clash_percentage",  default=10, help="maximum number of residues that clash between new loop and host, as a percentage of the length of the new loop (default is 10)" )
    parser.add_argument( "--loop_percentage",   default=50, help="maximum percentage of new loop as 'loop' secondary structure (default is 50)" )
    parser.add_argument( "--open_tolerances", help="set all tolerances to maximum", action="store_true" )
    # "--half-tolerances"; if passed, set all tolerances to half
    
    args = parser.parse_args()


    # maximize tolerances
    if args.open_tolerances:
        args.rmsd = 100
        args.length = 10000
        args.clash_percentage = 100
        args.loop_percentage = 100
        print 'tolerances have been set to maximum'

    # start pymol if appropriate
    if args.visualize:
        print 'importing pymol...'
        sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/')
        import pymol

        silence = pymol.finish_launching()
        silence = pymol.cmd.do('run ~/Desktop/PyRosetta/PyMOLPyRosettaServer.py')
        time.sleep(2)
    

    if True: #with open( 'thorp-log', 'w') as loggy, open('matching_domains.'+str(time.strftime("%Y%m%d-%H%M%S")), 'w') as matchy:

        orig_pdb_file = args.ref_pdb                # ./path/to/aaaa.pdb
        orig_pdb_id = orig_pdb_file.split('/')[-1]  # aaaa.pdb
        residue1 = args.ref_res_1
        residue2 = args.ref_res_2
        
        
        # DEBUG OUTPUT
        #output_string = '\n\n#######################################\n'+ \
        #    'This file contains debugging output for thorp.py\n'+ \
        #    'invoked at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'+ \
        #    'Reference PDB: '+str(orig_pdb_file)+'\n'+ \
        #    'Residues '+str(residue1)+' and '+str(residue2)+'\n' #+ \
        #    #str(len(os.listdir(family_dir)))+' dir items (PDBs to posify and traverse for transforms)\n\n'
        #loggy.write(output_string)
        #loggy.flush()

        # TRANSFORMS OUTPUT
        # TODO output tolerances
        #matchy_string = \
        #    'This list from *** directory, ### pdbs, *** tolerances, etc \n' + \
        #    'saved at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n'+ \
        #    str(sys.argv)+'\n\n\n'
        #matchy.write(matchy_string)
        #matchy.flush()
        
        find_matching_domains( args )

        #view_matching_domains()
        ## cycles through by default


