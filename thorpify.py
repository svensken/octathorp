from rosetta import *

init()

# reconstruct an RT object from a stored string
m = numeric.xyzMatrixdouble(1)
m = m.rows(11,12,13,14,15,16,17,18,19)
v = numeric.xyzVector_double(1,2,3)

rt = RT()
rt.set_rotation(m)
rt.set_translation(v)

#print rt
#print rt.get_translation()
#print rt.get_rotation()

all_s = ''
for life in os.listdir('alpha-beta-hydrolases/')[:10]:
    try:
        pose=Pose()
        pose_from_pdb(pose, 'alpha-beta-hydrolases/'+life)
    except:
        pass
    
    DsspMover().apply(pose)
    all_s +=  pose.secstruct()
