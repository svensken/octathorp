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


ap1_params = ['AP1.params']
type_set = generate_nonstandard_residue_set( ap1_params )
pose = pose_from_pdb( type_set, 'AP1_0001.pdb' )

tst_params = ['TST.params']
tst_type_set = generate_nonstandard_residue_set( tst_params )

tstpose = pose_from_pdb( tst_type_set, 'test.pdb')








