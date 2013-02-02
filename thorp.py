#!/usr/bin/env python

from rosetta import *
import os, sys, time, numpy, datetime

init()


# open file for log things
f = open('loggy', 'w')
raw_input("please run 'tail -f loggy' to see all script output, then press enter")



##################################################
#- Load all poses -------------------------------#
#TODO make inline functions
#TODO pickle the objects, though skipping the py part might take care of this
#TODO detect domains instead of manually specifying residues
## --> ExPASy's Prosite, Conserved Domain db, ESTHER db, ...
#TODO run entirely without input

try:
    pdb_file = [arg for arg in sys.argv if arg[-4:] == ".pdb"][0]
    # last two args are residues, hopefully
    residue1 = int(sys.argv[-2])
    residue2 = int(sys.argv[-1])
except:
    print "\nUsage: \n\
    $ thorp.py [file.pdb] [int residue 1] [int residue 2] \n\
    \n\
    Please put all relevant PDB files in the pdb directory.\n"
    #Please indicate path to PDB file relative to working directory
    sys.exit(0)

# load reference PDB
# (why does Pose error out?)
pose = Pose()
try:
    pose_from_pdb( pose, 'sample_pdbs/' + pdb_file )
    pose_load_result = 'no error'
except PyRosettaException:
    pose_load_result = 'slight error'

output_string = \
    'This file contains output from the thorp.py script\n'+ \
    'invoked at: '+str(time.strftime("%Y/%m/%d %H:%M:%S"))+'\n\n\n'+ \
    'Reference PDB: '+str(pdb_file)+'\n'+ \
    'Reference jump residues: '+str(residue1)+' and '+str(residue2)+'\n'+ \
    str(pose_load_result)+' loading the reference pose (no need to worry)\n'+ \
    'total reference residues: '+str(pose.total_residue())+'\n\n'
f.write(output_string)
f.flush()


# load all other PDB's
# cleaned pdb's end with 'A'
clean_pdb_file_list = [item for item in os.listdir('alpha-beta-hydrolases/') if item[-5:] == 'A.pdb']
output_string = '\n'+str(len(clean_pdb_file_list))+' PDBs to process into poses:\n'
f.write(output_string)
f.flush()

clean_pose_list = []
# load all the poses... O_o *gulp*...
n=0
t0 = time.time()
for thing in clean_pdb_file_list:
    n=n+1
    try:
        pose = Pose()
        pose_from_pdb( pose, 'alpha-beta-hydrolases/' + thing )
        pose_load_result = '.'
    except PyRosettaException:
        pose_load_result = '*'
    f.write(pose_load_result)
    if n%50==0 or n==len(clean_pdb_file_list):
        f.write(str(n))
    f.flush()
    # this is where we need to pickle the Pose objects
    #clean_pose_list.append(pose)
f.write( '\n\nAlles Klar\nTime taken to load all other poses: '+str(time.time()-t0) )
f.flush()
# writing to log files has never been uglier

#------------------------------------------------#
##################################################


"""

##################################################
#- Get all jumps --------------------------------#

# start with reference
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

pose_atom_tree = pose.atom_tree()
reference_transform = pose_atom_tree.get_stub_transform(stub1, stub2)
# break the numbers down
reference_jump = [ float(n) for n in str(reference_transform).split()[1:] ] # no 'RT'

# now all other jumps
new_pose = Pose()
try:
    pose_from_pdb(new_pose, 'sample_pdbs/' + pdb_file)
except PyRosettaException:
    print 'temp pose (', pdb_file, ') loaded'

stub_list = []
for middle_residue in range(2,new_pose.total_residue()):
    # first stub's central residue is #2
    # last stub's central residue is #last-1
    # Stub 1
    n_s1a1 = AtomID(1, middle_residue - 1)
    n_s1a2 = AtomID(1, middle_residue)
    n_s1a3 = AtomID(1, middle_residue + 1)
    n_stub = StubID(n_s1a1, n_s1a2, n_s1a3)
    stub_list.append(n_stub)
print 'stub_list built'

# build jump list
new_pose_atom_tree = new_pose.atom_tree()
jump_list = []
n = 0
for first_stub in stub_list:
    for second_stub in stub_list[ n: ]: # start at subsequent stub (to first_stub)
        new_jump_obj = new_pose_atom_tree.get_stub_transform(first_stub, second_stub)
        jump_list_element = [first_stub, second_stub, new_jump_obj]
        # structure of jump_list elements:
        # [ ..., [res1, res2, jump], ... ]
        jump_list.append(jump_list_element)
    n = n + 1
print 'jump_list built'

#------------------------------------------------#
##################################################



##################################################
#- Grade similarity of jumps to ref jump --------#

grade_list = []
t0 = time.time()
for jamp in jump_list:
    temp_jump = [ float(t) for t in str( jamp[2] ).split() [1:] ] # omit 'RT'
    ## Alternatives;
    # check if 0 first
    # (max-min)/max
    # shift all +10
    # (need a better way)
    big_diff = [ (a-b) / b for a,b in zip(temp_jump, reference_jump) ]
    avg_diff = numpy.average(big_diff)
    #print ''
    #print str(jamp[0]).split()[5]#10
    #print str(jamp[1]).split()[5]#10
    #print 'ref: ', reference_jump
    #print 'tmp: ', temp_jump
    #print 'diff list: ', big_diff
    #print 'avg diff:  ', avg_diff
    if abs(avg_diff) < .0005:
        grade_element = str(jamp[0]).split()[10], str(jamp[1]).split()[10], avg_diff
        #print grade_element
        grade_list.append(grade_element)

print '\ngrading time\n', time.time() - t0
sorted_grades = sorted(grade_list, key=lambda tup: abs(tup[2]))
print '\nSorted grades:'
for j in sorted_grades: # seperate lines
    print j

#------------------------------------------------#
##################################################



########
# deem loop's secondary composition
# how to confine? just limit string to residue numbers?
#ss = pose.secstruct()
#percent_helical = 100. * ss.count("H") / len(ss)
#H=helical, E=sheet, L=loop
#print str(percent_helical) + "% Helical"

####################################################
#### For Loop done #################################

#find_stub_transform (?)
# last but not least;
# http://graylab.jhu.edu/pyrosetta/downloads/scripts/demos/D120_Ligand_interface.py
# http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d7/d60/classnumeric_1_1xyz_vector.html
# http://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/de/d93/classnumeric_1_1xyz_matrix.html


"""

# close log file
f.close()
