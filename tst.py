from rosetta import *
import structural_alignment
reload(structural_alignment)

init()

pose1in=Pose()
pose2in=Pose()

pose_from_pdb( pose1in, 'alpha-beta-hydrolases/7yasA.pdb' )
pose_from_pdb( pose2in, 'alpha-beta-hydrolases/3puhA.pdb' )

r = range(33, 45)

structural_alignment.kabsch_alignment(pose1in, pose2in, r, r)
print '.'
