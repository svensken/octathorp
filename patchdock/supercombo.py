#!/usr/bin/env python

import os, math


def make_my_params( bot, cap, dist_line1, dist_line2 ):
    return """
########################################################################
#    PatchDock Parameter File for """+ bot +'-'+ cap +"""
########################################################################

#    File Names:
receptorPdb """+bot+"""
ligandPdb """+cap+"""
receptorMs """+bot+""".ms
ligandMs """+cap+""".ms
#receptorActiveSite site1.txt
#ligandActiveSite site2.txt
protLib /u1/home/svensken/PatchDock/chem.lib
log-file patch_dock.log
log-level 2

#   Distance constraint parameters:
#       distanceConstraints <rec_atom_index> <lig_atom_index> <dist_thr>
# <rec_atom_index> - receptor atom used for constraint
# <lig_atom_index> - ligand atom used for constraint
# <dist_thr> - maximum allowed distance between the specified atom centers
#distanceConstraints rec_atom_index lig_atom_index dist_thr
"""             + \
dist_line1      + \
'\n'            + \
dist_line2      + \
"""
#distanceConstraintsFile dists.cst

#    Surface Segmentation Parameters:
#       receptorSeg <low_patch_thr> <high_patch_thr> <prune_thr>
#                   <knob> <flat> <hole>
#                   <hot spot filter type>
#    <low_patch_thr>,<high_patch_thr> - min and max patch diameter
#    <prune_thr> - minimal distance between points inside the patch
#    <knob> <flat> <hole> - types of patches to dock (1-use, 0-do not use) (may need tuning)
#    <hot spot filter type> :None - 0, Antibody - 1, Antigen - 2
#                             Protease - 3, Inhibitor - 4, Drug - 5
receptorSeg 10.0 20.0 1.5 1 0 1 0
ligandSeg 10.0 20.0 1.5 1 0 1 0

#    Scoring Parameters:
#        scoreParams <small_interfaces_ratio> <max_penetration> <ns_thr>
#                    <rec_as_thr> <lig_as_thr> <patch_res_num> <w1 w2 w3 w4 w5>
#    <small_interfaces_ratio> - the ratio of the low scoring transforms to be removed
#    <max_penetration> - maximal allowed penetration between molecules surfaces
#    <ns_thr> - normal score threshold
#    <rec_as_thr> <lig_as_thr> - the minimal ratio of the active site area in the solutions
#    <patch_res_num> - number of results to consider in each patch
#    <w1 w2 w3 w4 w5> - scoring weights for ranges:
#                [-5.0,-3.6],[-3.6,-2.2],[-2.2,-1.0],[-1.0,1.0],[1.0-up]
scoreParams 0.3 -5.0 0.5 0.0 0.0 1500 -8 -4 0 1 0

#    Desolvation Scoring Parameters:
#        desolvationParams <energy_thr> <cut_off_ratio>
#    <energy_thr> - remove all results with desolvation energy higher than threshold
#    <cut_off_ratio> - the ratio of low energy results to be kept
#    First filtering with energy_thr is applied and the remaining results
#    can be further filtered with cut_off_ratio.
desolvationParams 500.0 1.0

########################################################################
#   Advanced Parameters
########################################################################

#    Clustering Parameters:
#    clusterParams < rotationVoxelSize > < discardClustersSmaller > < rmsd > < final clustering rmsd >
clusterParams 0.1 4 2.0 4.0

#    Base Parameters:
#    baseParams <min_base_dist> <max_base_dist> <# of patches for base: 1 or 2>
baseParams 4.0 13.0 2

#    Matching Parameters:
#  matchingParams <geo_dist_thr> <dist_thr> <angle_thr> <torsion_thr> 
#     <angle_sum_thr>
matchingParams 1.5 1.5 0.4 0.5 0.9
# 1 - PoseClustering (default), 2 - Geometring Hashing
matchAlgorithm 1

#    Grid Parameters:
#      receptorGrid <gridStep> <maxDistInDistFunction> <vol_func_radius>
#      NOTE: the vol_func_radius of small molecules and peptides should be 3A!
receptorGrid 0.5 6.0 6.0
ligandGrid 0.5 6.0 6.0

# Energy Parameters:
vdWTermType 1
attrVdWEnergyTermWeight 1.01
repVdWEnergyTermWeight 0.5
HBEnergyTermWeight 1.0
ACE_EnergyTermWeight 1.0
piStackEnergyTermWeight 0.0
confProbEnergyTermWeight 0.1
COM_distanceTermWeight 1.07
energyDistCutoff 6.0
elecEnergyTermWeight 0.1
radiiScaling 0.8
"""


def distance(x1, y1, z1, x2, y2, z2):
    dist = math.sqrt( (float(x1)-float(x2))**2 + (float(y1)-float(y2))**2 + (float(z1)-float(z2))**2 )
    return dist



#path
p = '/home/svensken/octathorp/patchdock/'

bottoms = {'1WM1bot.pdb':'138_248', '2WUFbot.pdb':'139_222', '3ANSbot.pdb':'361_485', '3B12bot.pdb':'128_230', '3GZJbot.pdb':'111_198'}
caps = {'1WM1cap.pdb':'140_238', '2WUFcap.pdb':'145_216', '3ANScap.pdb':'369_474', '3B12cap.pdb':'139_226', '3GZJcap.pdb':'121_195'}

#distanceConstraints 1022 1075 16.3
#distanceConstraints 1030 1783 16.9
#ATOM   1058  CA  ILE A 138      20.905  -0.983   1.271  1.00 24.12           C


for bot in bottoms:
    for cap in caps:

        os.mkdir( 'combos/'+bot+'_'+cap )
        os.chdir( 'combos/'+bot+'_'+cap )

        with open( p+bot, 'r' ) as mr_bottom:
            lines = mr_bottom.readlines()
            res_nums = bottoms[bot].split('_')
            res1 = res_nums[0]
            res2 = res_nums[1]
            for line in lines:
                try:
                    if line.split()[5]==res1 and line.split()[2]=='CA':
                        bot_atomid1 = line.split()[1]
                        bot_x1 = line.split()[6]
                        bot_y1 = line.split()[7]
                        bot_z1 = line.split()[8]
                    if line.split()[5]==res2 and line.split()[2]=='CA':
                        bot_atomid2 = line.split()[1]
                        bot_x2 = line.split()[6]
                        bot_y2 = line.split()[7]
                        bot_z2 = line.split()[8]
                except:
                    pass

        with open( p+cap, 'r' ) as mr_cap:
            lines = mr_cap.readlines()
            res_nums = caps[cap].split('_')
            res1 = res_nums[0]
            res2 = res_nums[1]
            for line in lines:
                try:
                    if line.split()[5]==res1 and line.split()[2]=='CA':
                        cap_atomid1 = line.split()[1]
                        cap_x1 = line.split()[6]
                        cap_y1 = line.split()[7]
                        cap_z1 = line.split()[8]
                    if line.split()[5]==res2 and line.split()[2]=='CA':
                        cap_atomid2 = line.split()[1]
                        cap_x2 = line.split()[6]
                        cap_y2 = line.split()[7]
                        cap_z2 = line.split()[8]
                except:
                    pass


        dist1 = distance( bot_x1, bot_y1, bot_z1, 
                          cap_x1, cap_y1, cap_z1 )
        dist2 = distance( bot_x2, bot_y2, bot_z2, 
                          cap_x2, cap_y2, cap_z2 )
        
        cst_line1 = 'distanceConstraints '+bot_atomid1+' '+cap_atomid1+' '+str(dist1)
        cst_line2 = 'distanceConstraints '+bot_atomid2+' '+cap_atomid2+' '+str(dist2)

        with open('params.txt','w') as params:
            params_template = make_my_params( p+bot, p+cap, cst_line1, cst_line2 )
            params.write( params_template )

        # is this specific to this pair or can this be run on individual pdbs beforehand?
        os.system('buildMS.pl '+ p+bot +' '+ p+cap ) 

        os.system('patch_dock.Linux params.txt output.txt' ) 

        os.system('transOutput.pl output.txt 1 10000') # output all

        os.chdir( p )

        print 'klar med en'
