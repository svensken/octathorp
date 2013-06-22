#!/usr/bin/env python

from rosetta import *
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res']
args = utility.vector1_string()
args.extend( opts )
core.init( args )



pose = Pose()
pose_from_pdb( pose, '3ANS_trimmed.pdb' )


pi = pose.pdb_info()
print 'pose'
for i in [361, 369, 474, 485]: 
    print i, pi.pdb2pose('A', i)
    print i, pi.pdb2pose('B', i)









