#!/usr/bin/env python

import os
from rosetta import *
from rosetta.protocols.rigid import *
import structural_alignment

opts = ['app', 
        '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), 
        '-ignore_unrecognized_res', 
        '-extra_res_fa','/home/svensken/octathorp/deepak/AP1.params', 
                        '/home/svensken/octathorp/deepak/AP2.params',
                        '/home/svensken/octathorp/deepak/AP3.params',
                        '/home/svensken/octathorp/deepak/TST.params' ]

args = utility.vector1_string()
args.extend( opts )
core.init( args )




# get ligand coordinates
ap1 = Pose()
pose_from_pdb(ap1, 'AP1_0001.pdb')
ap1_c33 = ap1.residue(1).atom('C33').xyz()
ap1_s1  = ap1.residue(1).atom('S1' ).xyz()
ap1_o9  = ap1.residue(1).atom('O9' ).xyz()
ap1_c29 = ap1.residue(1).atom('C29').xyz()

ap2 = Pose()
pose_from_pdb(ap2, 'AP2_0001.pdb')
ap2_c1  = ap2.residue(1).atom('C1' ).xyz()
ap2_s1  = ap2.residue(1).atom('S1' ).xyz()
ap2_o1  = ap2.residue(1).atom('O1' ).xyz()
ap2_c3  = ap2.residue(1).atom('C3' ).xyz()

ap3 = Pose()
pose_from_pdb(ap3, 'AP3_0001.pdb')
ap3_c1  = ap3.residue(1).atom('C1' ).xyz()
ap3_s1  = ap3.residue(1).atom('S1' ).xyz()
ap3_o1  = ap3.residue(1).atom('O1' ).xyz()
ap3_c3  = ap3.residue(1).atom('C3' ).xyz()


for filename in os.listdir('/home/deegupta/matching/matches_all'):
    filename = '/home/deegupta/matching/matches_all/'+filename

    with open(filename, 'r') as pdbfile:
        lines = pdbfile.readlines()
        for line in lines:
            terms = line.split()
            
            try:
                if terms[2]=='C1' and terms[3]=='TST' and terms[4]=='X':
                    tst_c1 = [ terms[6], terms[7], terms[8] ]
                if terms[2]=='S1' and terms[3]=='TST' and terms[4]=='X':
                    tst_s1 = [ terms[6], terms[7], terms[8] ]
                if terms[2]=='O1' and terms[3]=='TST' and terms[4]=='X':
                    tst_o1 = [ terms[6], terms[7], terms[8] ]
                if terms[2]=='C2' and terms[3]=='TST' and terms[4]=='X':
                    tst_c2 = [ terms[6], terms[7], terms[8] ]
                    break
            except:
                pass

    pose = Pose()
    pose_from_pdb( pose, filename )

    ap1, pose = structural_alignment.kabsch_alignment(  ap1, 
                                                        pose, 
                                                        [   ap1_c33,
                                                            ap1_s1,
                                                            ap1_o9,
                                                            ap1_c29  ], 
                                                        [   tst_c1,
                                                            tst_s1,
                                                            tst_o1,
                                                            tst_c2   ] )


    print ap1
    print pose

    #dump_pdb pose
    #dump_pdb pep

    #system cat pepfile >> posefile

