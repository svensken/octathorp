#!/usr/bin/env python

import time, sys, os

print 'importing rosetta...'
from rosetta import *
from rosetta.protocols.rigid import *
# silence output
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res', '-extra_res_cen', '/home/svensken/octathorp/molfiles/1WM1_PTB.cen.params', '-extra_res_fa', '/home/svensken/octathorp/molfiles/1WM1_PTB.fa.params' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )

#print 'importing pymol...'
#sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/2.6/lib/python2.6/site-packages/')
# silence output, no GUI
#import __main__
#__main__.pymol_argv = ['pymol','-qc']
#import pymol
#pymol.finish_launching()
#pymol.cmd.do('run ~/Desktop/PyRosetta/PyMOLPyRosettaServer.py')

this_dir = os.path.basename( os.getcwd() )

def ce_align_poses():
    pymover = PyMOL_Mover()
    pymover.apply(hpose)
    pymover.apply(gpose)
    #time.sleep(1)
    pymol.cmd.load( hpose_filename )
    pymol.cmd.load( gpose_filename )
    hpose_name = '4fwb'
    gpose_name = '2yas'
    pymol.cmd.cealign(hpose_name, gpose_name)
    pymol.cmd.save('ce_host_temp.pdb', hpose_name)
    pymol.cmd.save('ce_guest_temp.pdb', gpose_name)
    #pose = Pose()
    #pose = Pose()
    pose_from_pdb(hpose, 'ce_host_temp.pdb')
    pose_from_pdb(gpose, 'ce_guest_temp.pdb')
    pymover.apply(hpose)
    pymover.apply(gpose)

    return hpose, gpose

def setup_ligand_params(): 
    # no workie
    params_list = Vector1( ['molfiles/1WM1_PTB.cen.params','molfiles/1WM1_PTB.fa.params'] )

    res_set = ChemicalManager.get_instance().nonconst_residue_type_set('fa_standard')
    atoms = ChemicalManager.get_instance().atom_type_set('fa_standard')
    mm_atoms = ChemicalManager.get_instance().mm_atom_type_set('fa_standard')
    orbitals = ChemicalManager.get_instance().orbital_type_set('fa_standard')
    elements = ChemicalManager.get_instance().element_set('fa_standard')

    res_set.read_files(params_list, atoms, elements, mm_atoms, orbitals)

    return res_set

def find_termini(pose, hres1):
    pi = pose.pdb_info()
    first_hit_B = True
    r = 1
    while r <= pose.total_residue():    
        # detect ends of chain B
        if pi.chain(r) == 'B':
            if first_hit_B:
                first_B = r
                first_hit_B = False
            last_B = r
        r += 1

    hterm1 = hres1
    hterm2 = hres1+1
    gterm1 = first_B
    gterm2 = last_B

    return hterm1, hterm2, gterm1, gterm2

def snip_poses(hpose,hres1,hres2,gpose,gres1,gres2):
    # delete between hres1 hres2
    # list should be [ hres2-1, hres2-2, ... , hres1+1 ]
    for r in reversed(range(hres1+1, hres2)): # WILL SEGFAULT AT 3r = 1
        hpose.delete_polymer_residue( r )

    # add guest domain to hpose as new chain
    #append_residue_by_jump(hpose, gpose.res, hres1, start_new_chain=True)
    hpose.append_residue_by_jump(gpose.residue(gres1), hres1, start_new_chain=True)
    pi=hpose.pdb_info()
    for r in reversed(range(gres1, gres2+1)):
        hpose.append_polymer_residue_after_seqpos(gpose.residue(r), hpose.total_residue(), 0)
        pi.chain( hpose.total_residue(), "B" )

    # reload to make new chain behave correctly
    hpose.dump_pdb('reload_temp.pdb')
    pose=Pose()
    pose_from_pdb(pose, 'reload_temp.pdb')

    return pose

def setup_pose():
    if pdbs_prepared:
        pose = Pose()
        pose_from_pdb( pose, sys.argv[1] )
        # for benchmark, let's just define them manually
        # POSE NUMBERING
        hterm1 = 137 #140
        hterm2 = 138 #251
        gterm1 = 205
        gterm2 = 310

        #if nonstandard_ligand:
        #    res_set = setup_ligand_params()
        #    # sys.argv = script.py mixture.pdb
        #    pose_from_pdb( pose, res_set, sys.argv[1] )

    else:
        # sys.argv = script.py host hcut1 hcut2 guest gcut1 gcut2
        hpose_filename = sys.argv[1]
        gpose_filename = sys.argv[4]
        hpose = Pose()
        pose_from_pdb(hpose, hpose_filename)
        gpose = Pose()
        pose_from_pdb(gpose, gpose_filename)
        # convert to pose numbering
        hpi = hpose.pdb_info()
        gpi = gpose.pdb_info()
        # sys.argv = script.py host hcut1 hcut2 guest gcut1 gcut2
        hres1 = int( sys.argv[2] )
        hres2 = int( sys.argv[3] )
        gres1 = int( sys.argv[5] )
        gres2 = int( sys.argv[6] )
        hres1 = hpi.pdb2pose( "A", hres1 )
        hres2 = hpi.pdb2pose( "A", hres2 )
        gres1 = gpi.pdb2pose( "A", gres1 )
        gres2 = gpi.pdb2pose( "A", gres2 )

        #hpose, gpose = ce_align_poses()
        
        pose = snip_poses(hpose,hres1,hres2,gpose,gres1,gres2)

        hterm1, hterm2, gterm1, gterm2 = find_termini(pose, hres1)

    hterm1 = AtomID( 1, hterm1 ) #hres1
    hterm2 = AtomID( 1, hterm2 ) #hres2
    gterm1 = AtomID( 1, gterm1 )#int(str(pose.fold_tree().jump_edge(1)).split()[2]) ) #chainB first res
    gterm2 = AtomID( 1, gterm2 )#pose.total_residue()) #chainB last res

    return hterm1, hterm2, gterm1, gterm2, pose

def setup_movers_and_sf():
    pass

def setup_constraints():
    pass
    

##################################################


pdbs_prepared = True
#nonstandard_ligand = False # works better by passing extra_res flags to init()
jump_number = 2

hterm1, hterm2, gterm1, gterm2, pose = setup_pose()

print hterm1
print hterm2
print gterm1
print gterm2

#sidechain_pose = Pose()
#sidechain_pose.assign(pose)
#to_centroid = SwitchResidueTypeSetMover('centroid')
#to_centroid.apply(pose)
to_centroid = SwitchResidueTypeSetMover('centroid')
to_fullatom = SwitchResidueTypeSetMover('fa_standard')

recover_sidechains = ReturnSidechainMover(pose)

to_centroid.apply(pose)
centroid_pose = Pose()
centroid_pose.assign(pose)

# SETUP MOVERS & SCOREFXN
randomize_upstream = RigidBodyRandomizeMover(pose, jump_number, partner_upstream)
randomize_downstream = RigidBodyRandomizeMover(pose, jump_number, partner_downstream)
pert_mover = RigidBodyPerturbMover( jump_number, 3, 8 ) #(jump_num, translation, rotation)
spin = RigidBodySpinMover( jump_number )
slide_into_contact = DockingSlideIntoContact( jump_number )

#sf = create_score_function('interchain_cen')
scorefxn_low = create_score_function('interchain_cen')
scorefxn_low.set_weight( atom_pair_constraint, 5 )
scorefxn_high = create_score_function('docking')
scorefxn_high.set_weight( atom_pair_constraint, 5 )
scorefxn_high_min = create_score_function_ws_patch('docking', 'docking_min')
scorefxn_high_min.set_weight( atom_pair_constraint, 5 )

movemap = MoveMap()
movemap.set_jump( jump_number, True )
minmover = MinMover()
minmover.movemap(movemap)
minmover.score_function(scorefxn_high)

perturb = SequenceMover()
perturb.add_mover(randomize_upstream)
perturb.add_mover(randomize_downstream)
perturb.add_mover(pert_mover)
perturb.add_mover(spin)
#perturb.add_mover(slide_into_contact)
perturb.add_mover(to_fullatom)
perturb.add_mover(recover_sidechains)
perturb.add_mover(minmover)
#


# CONSTRAINTS
GF1 = constraints.GaussianFunc( 14.1, 2.0 )
GF2 = constraints.GaussianFunc( 11.5, 2.0 )
apc1 = constraints.AtomPairConstraint( hterm1, gterm1, GF1 )
apc2 = constraints.AtomPairConstraint( hterm2, gterm2, GF2 )
pose.add_constraint( apc1 )
pose.add_constraint( apc2 )



#scorefxn_high.show(pose)


# DOCKER & JOB-DISTRIBUTOR

dock_prot = DockingProtocol()    # contains many docking functions
dock_prot.set_movable_jumps(Vector1([ jump_number ]))
dock_prot.set_lowres_scorefxn(scorefxn_low)
dock_prot.set_highres_scorefxn(scorefxn_high_min)


#AddPyMolObserver(pose, True)


# this not needed if we start with a clean 30 directories every time
#if os.path.exists('energies.csv'):
#    with open('energies.sc', 'w') as csv:
#        csv.write( '' )


#jd = PyJobDistributor('jd_output', 400, sf)
#jd.native_pose = pose

#while not jd.job_complete:
for i in range(100):
    
    filename = 'manual_' + str(i) + '.pdb'
    print filename

    centroid_pose.assign(pose)

    # perturb the structure
    perturb.apply(centroid_pose)

    print 'perturbed'

    # perform docking
    #sf.show(pose)
    dock_prot.apply(centroid_pose)
    #sf.show(pose)

    print 'docked'

    # record scores
    cen_score_with_cst = str(scorefxn_low.score(centroid_pose))
    scorefxn_low.set_weight( atom_pair_constraint, 0 )
    cen_score_without_cst = str(scorefxn_low.score(centroid_pose))

    print 'scores saved'

    # write scores
    with open('energies.sc','a') as scorefile:
        linetowrite = cen_score_with_cst+','+cen_score_without_cst+','+this_dir+'/'+filename
        scorefile.write( linetowrite + '\n' )

    print 'written'

    #recover_sidechains.apply(pose)
    to_fullatom.apply(centroid_pose) 

    print 'fullatomed'

    #jd.output_decoy(pose)

    # dump fullatom pdb (manually)
    centroid_pose.dump_pdb( filename ) #dump_scored_pdb( filename, sf )

    print filename, 'dumped'

    with open('status.update', 'a') as statusupdate:
        statusupdate.write(str(time.strftime("%Y/%m/%d-%H:%M:%S"))+'\n')
