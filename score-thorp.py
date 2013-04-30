#!/usr/bin/env python

print "importing rosetta"
from rosetta import *
import os

opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )

# parameters
hpose = Pose()
pose_from_pdb(hpose, './4fwb.pdb')
gpose = Pose()
pose_from_pdb(gpose, './2yas.pdb')
hres1 = 136 #remains 136
hres2 = 212 #becomes 137
gres1 = 113 #becomes 219
gres2 = 185 #becomes 291


# rmsd ref
rmsd_ref_pose = Pose()
pose_from_pdb( rmsd_ref_pose, './bashy/42/output_1.pdb' ) #lowest score



# run from octathorp/ and search through bashy/##'s
for root, dirnames, filenames in os.walk('./bashy'):
    for f in filenames:
        if f.startswith('output_') and f.endswith('.pdb'):
            this_pdb = os.path.join( root, f )
            
            pose = Pose()
            pose_from_pdb( pose, this_pdb )

            pi = pose.pdb_info()

            # assume chB goes res r to total_residue()
            r = 1
            while pi.chain(r) != 'B':
                r += 1
            hterm1_res = pose.residue( hres1 )   #136
            hterm2_res = pose.residue( hres1+1 ) #137
            
            gterm1_res = pose.residue( r+1 )
            gterm2_res = pose.residue( pose.total_residue() )

            # DISTANCES
            dist1 = hterm1_res.xyz('CA').distance( gterm1_res.xyz('CA') )
            dist2 = hterm2_res.xyz('CA').distance( gterm2_res.xyz('CA') )

            print 
            print this_pdb
            print hterm1_res.seqpos(), hterm2_res.seqpos(), gterm1_res.seqpos(), gterm2_res.seqpos()
            print dist1, dist2

            # TOTAL ENERGY
            # fullatom
            fa_e = pose.energies()
            fa_total_energy = fa_e.total_energy()

            # centroid
            to_centroid = SwitchResidueTypeSetMover('centroid')
            to_centroid.apply(pose)

            c_e = pose.energies()
            c_total_energy = c_e.total_energy()

            #print 'fa',fa_total_energy
            #print 'c ', c_total_energy


            sf = create_score_function('interchain_cen')
            sf.set_weight( atom_pair_constraint, 10 )

            GF = constraints.HarmonicFunc( 6.0, 2.0 ) #GaussianFunc( 6.0, 2.0 )
            apc1 = constraints.AtomPairConstraint( AtomID(1,hres1), AtomID(1,r + 1), GF ) # hterm1, gterm1
            apc2 = constraints.AtomPairConstraint( AtomID(1,hres1+1), AtomID(1,pose.total_residue()), GF ) # hterm2, gterm2
            pose.add_constraint( apc1 )
            pose.add_constraint( apc2 )

            # with our constraints
            e = pose.energies()
            te = e.total_energy()
            #scores_list.append(te)

            # default scoring
            sf.set_weight( atom_pair_constraint, 0 )
            d_e = pose.energies()
            d_te = e.total_energy()


            # RMSD
            rmsd = CA_rmsd( rmsd_ref_pose, pose, r+1, pose.total_residue() )

            # record to csv
            # open and close every loop for a write every time
            with open('rmsd_vs_energy.csv', 'a') as csv:
                csv.write( str(rmsd)+', '+str(te)+', '+str(d_te)+'\n' )
