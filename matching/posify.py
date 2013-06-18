#!/usr/bin/env python

from rosetta import *
opts = [ 'app', '-database', os.path.abspath( os.environ['PYROSETTA_DATABASE'] ), '-mute', 'all', '-ignore_unrecognized_res']#, '-extra_res_cen', '/home/svensken/octathorp/molfiles/1WM1_PTB.cen.params', '-extra_res_fa', '/home/svensken/octathorp/matching/params/EVS.fa.params', '/home/svensken/octathorp/matching/params/S38.fa.params', '/home/svensken/octathorp/matching/params/FAH.fa.params', '/home/svensken/octathorp/matching/params/old.KEM.fa.params', '/home/svensken/octathorp/matching/params/KEM.fa.params', '/home/svensken/octathorp/matching/params/PTB.fa.params' ]
args = utility.vector1_string()
args.extend( opts )
core.init( args )

#wm = Pose()
#pose_from_pdb( wm, 'new.1WM1.pdb' )

w1 = Pose()
pose_from_pdb( w1, 'natives/3GZJ.pdb')#'natives/new.2WUF.pdb' )


#wu = Pose()
#pose_from_pdb( wu, 'natives/new.2WUF.pdb' )

#an = Pose()
#pose_from_pdb( an, 'natives/new.3ANS.pdb' )
#
#b1 = Pose()
#pose_from_pdb( b1, 'natives/new.3B12.pdb' )
#
#gz = Pose()
#pose_from_pdb( gz, 'natives/3GZJ.pdb' )



#pi = wm.pdb_info()
#print 'wm'
#for i in [1,2,3,4]:
#    print i, pi.pdb2pose('A', i)

p1 = w1.pdb_info()
print 'wu_an, wu_b1'
for i in [87, 244, 216, 128]: # 269 too far away
    print i, p1.pdb2pose('A', i)
#pi = wu.pdb_info()
#print 'wu'
#for i in [192, 54, 270, 113]: # 269 too far away
#    print i, pi.pdb2pose('A', i)
#
#pi = an.pdb_info()
#print 'an', an.total_residue()
#for i in [335, 524, 496, 383, 466]:
#    print i, pi.pdb2pose('A', i)
#
#pi = b1.pdb_info()
#print 'b1'
#for i in [104, 271, 128, 105, 108, 150, 149, 212]:
#    print i, pi.pdb2pose('A', i)
#
#pi = gz.pdb_info()
#print 'gz'
#for i in [87, 244, 216, 128, 125]:
#    print i, pi.pdb2pose('A', i)











