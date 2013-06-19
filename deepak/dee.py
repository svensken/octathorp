#!/usr/bin/env python

import os
from rosetta import *
from rosetta.protocols.rigid import *
import structural_alignment

opts = ['app', 
        '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), 
        '-ignore_unrecognized_res' ]

args = utility.vector1_string()
args.extend( opts )
core.init( args )


tst_type_set = generate_nonstandard_residue_set( ['TST.params'] )


# get ligand coordinates
ap1_type_set = generate_nonstandard_residue_set( ['AP1.params'] )
ap1 = pose_from_pdb(ap1_type_set, 'AP1_0001.pdb')
ap1_c33 = ap1.residue(1).atom('C33').xyz()
ap1_s1  = ap1.residue(1).atom('S1' ).xyz()
ap1_o9  = ap1.residue(1).atom('O9' ).xyz()
ap1_c29 = ap1.residue(1).atom('C29').xyz()
# listify
ap1_c33 = [ ap1_c33[0], ap1_c33[1], ap1_c33[2] ]
ap1_s1  = [ ap1_s1[0], ap1_s1[1], ap1_s1[2] ]
ap1_o9  = [ ap1_o9[0], ap1_o9[1], ap1_o9[2] ]
ap1_c29 = [ ap1_c29[0], ap1_c29[1], ap1_c29[2] ]

ap2_type_set = generate_nonstandard_residue_set( ['AP2.params'] )
ap2 = pose_from_pdb(ap2_type_set, 'AP2_0001.pdb')
ap2_c1  = ap2.residue(1).atom('C1' ).xyz()
ap2_s1  = ap2.residue(1).atom('S1' ).xyz()
ap2_o1  = ap2.residue(1).atom('O1' ).xyz()
ap2_c3  = ap2.residue(1).atom('C3' ).xyz()
# listify
ap2_c1  = [ ap2_c1[0], ap2_c1[1], ap2_c1[2] ]
ap2_s1  = [ ap2_s1[0], ap2_s1[1], ap2_s1[2] ]
ap2_o1  = [ ap2_o1[0], ap2_o1[1], ap2_o1[2] ]
ap2_c3  = [ ap2_c3[0], ap2_c3[1], ap2_c3[2] ]

ap3_type_set = generate_nonstandard_residue_set( ['AP3.params'] )
ap3 = pose_from_pdb(ap3_type_set, 'AP3_0001.pdb')
ap3_c1  = ap3.residue(1).atom('C1' ).xyz()
ap3_s1  = ap3.residue(1).atom('S1' ).xyz()
ap3_o1  = ap3.residue(1).atom('O1' ).xyz()
ap3_c3  = ap3.residue(1).atom('C3' ).xyz()
# listify
ap3_c1  = [ ap3_c1[0], ap3_c1[1], ap3_c1[2] ]
ap3_s1  = [ ap3_s1[0], ap3_s1[1], ap3_s1[2] ]
ap3_o1  = [ ap3_o1[0], ap3_o1[1], ap3_o1[2] ]
ap3_c3  = [ ap3_c3[0], ap3_c3[1], ap3_c3[2] ]



for filename in os.listdir('/home/deegupta/matching/matches_all'):
    # don't recycle!!!!!!
    ap1_tmp = None
    ap1_tmp = Pose()
    ap1_tmp.assign(ap1)
    ap2_tmp = None
    ap2_tmp = Pose()
    ap2_tmp.assign(ap2)
    ap3_tmp = None
    ap3_tmp = Pose()
    ap3_tmp.assign(ap3)

    pathname = '/home/deegupta/matching/matches_all/'+filename

    with open(pathname, 'r') as pdbfile:
        lines = pdbfile.readlines()
        for line in lines:
            terms = line.split()
            
            try:
                if terms[2]=='C1' and terms[3]=='TST' and terms[4]=='X':
                    tst_c1 = [ float(terms[6]), float(terms[7]), float(terms[8]) ]
                if terms[2]=='S1' and terms[3]=='TST' and terms[4]=='X':
                    tst_s1 = [ float(terms[6]), float(terms[7]), float(terms[8]) ]
                if terms[2]=='O1' and terms[3]=='TST' and terms[4]=='X':
                    tst_o1 = [ float(terms[6]), float(terms[7]), float(terms[8]) ]
                if terms[2]=='C2' and terms[3]=='TST' and terms[4]=='X':
                    tst_c2 = [ float(terms[6]), float(terms[7]), float(terms[8]) ]
                    break
            except:
                pass
    
    pose = pose_from_pdb( tst_type_set, pathname )
    

    ap1_tmp, pose1 = structural_alignment.kabsch_alignment( ap1_tmp, 
                                                        pose, 
                                                        [   ap1_c33,
                                                            ap1_s1,
                                                            ap1_o9,
                                                            ap1_c29  ], 
                                                        [   tst_c1,
                                                            tst_s1,
                                                            tst_o1,
                                                            tst_c2   ] )
    ap1_out = None
    ap1_out = Pose()
    ap1_out.assign(ap1_tmp)
    pose = pose_from_pdb( tst_type_set, pathname )
    ap2_tmp, pose2 = structural_alignment.kabsch_alignment( ap2_tmp, 
                                                        pose, 
                                                        [   ap2_c1,
                                                            ap2_s1,
                                                            ap2_o1,
                                                            ap2_c3  ], 
                                                        [   tst_c1,
                                                            tst_s1,
                                                            tst_o1,
                                                            tst_c2   ] )
    ap2_out = None
    ap2_out = Pose()
    ap2_out.assign(ap2_tmp)
    pose = pose_from_pdb( tst_type_set, pathname )
    ap3_tmp, pose3 = structural_alignment.kabsch_alignment( ap3_tmp, 
                                                        pose, 
                                                        [   ap3_c1,
                                                            ap3_s1,
                                                            ap3_o1,
                                                            ap3_c3  ], 
                                                        [   tst_c1,
                                                            tst_s1,
                                                            tst_o1,
                                                            tst_c2   ] )
    ap3_out = None
    ap3_out = Pose()
    ap3_out.assign(ap3_tmp)
    pose = pose_from_pdb( tst_type_set, pathname )

    #print ap1_tmp
    #print pose

    dump_pdb(ap1_out, 'tmp/ap1.pdb')
    dump_pdb(ap2_out, 'tmp/ap2.pdb')
    dump_pdb(ap3_out, 'tmp/ap3.pdb')
    dump_pdb(pose1, 'tmp/pose1.pdb')
    dump_pdb(pose2, 'tmp/pose2.pdb')
    dump_pdb(pose3, 'tmp/pose3.pdb')


    os.system('cat tmp/pose1.pdb tmp/ap1.pdb >> results/'+filename+'.ap1.pdb')
    os.system('cat tmp/pose2.pdb tmp/ap2.pdb >> results/'+filename+'.ap2.pdb')
    os.system('cat tmp/pose3.pdb tmp/ap3.pdb >> results/'+filename+'.ap3.pdb')

