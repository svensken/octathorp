#!/usr/bin/env python

import os, sys, subprocess


temp_resnum_dict = { '1WM1' : { 'botA' : 136,#138,
                                'capA' : 142,#140,
                                'capB' : 240,#238,
                                'botB' : 246},#248 },
                     '2WUF' : { 'botA' : 139,#139,
                                'capA' : 145,#145,
                                'capB' : 216,#216,
                                'botB' : 222},#222 },
                     '3ANS' : { 'botA' : 362,#361,
                                'capA' : 368,#369,
                                'capB' : 477,#474,
                                'botB' : 483},#485 },
                     '3B12' : { 'botA' : 131,#128,
                                'capA' : 137,#139,
                                'capB' : 225,#226,
                                'botB' : 231},#230 },
                     '3GZJ' : { 'botA' : 113,#111,
                                'capA' : 119,#121,
                                'capB' : 194,#195,
                                'botB' : 200}}#198 } }

temp_pose_num_dict = {}


### REORDER
os.chdir('matched/')
for match in os.listdir('.'):
  print 'will reorder: '+match
  pose_num = 1
  this_id = match[:4]
  # need to start with empty dict
  temp_pose_num_dict[this_id] = {}
  
  straight_As = []
  chunk_o_Bs = []
  with open(match, 'r') as pdb_in:
    lines = pdb_in.readlines()
    for line in lines:
      linelist = list(line)
      if line[:4] == "ATOM":
        res_num = int(line[23:26])
        if linelist[21]=='A' and line[12:15]==' N ' and res_num==temp_resnum_dict[this_id]['botB']:
          #print 'prev: ',prev_resi
          #print 'resi: ',resi
          #print 'delta: ', delta
          straight_As.append('repwace_me')
          straight_As.append(''.join(linelist))
        elif linelist[21] == 'B':
          linelist[21] = "A"
          chunk_o_Bs.append(''.join(linelist))
        else:
          straight_As.append(''.join(linelist))
      

      if line[:6] == "ENDMDL":
        break
  
  straight_As = ''.join(straight_As)
  chunk_o_Bs  = ''.join(chunk_o_Bs)
  straight_As = straight_As.replace('repwace_me', chunk_o_Bs)
  
  # store pose numbers
  for line in straight_As.split('\n'):
    if line[:4] == "ATOM" and line[12:15]==' N ':
      res_num = int(line[23:26])
      if res_num == temp_resnum_dict[this_id]['botA']:
        temp_pose_num_dict[this_id]['botA'] = pose_num
      if res_num == temp_resnum_dict[this_id]['capA']:
        temp_pose_num_dict[this_id]['capA'] = pose_num
      if res_num == temp_resnum_dict[this_id]['capB']:
        temp_pose_num_dict[this_id]['capB'] = pose_num
      if res_num == temp_resnum_dict[this_id]['botB']:
        temp_pose_num_dict[this_id]['botB'] = pose_num
      pose_num += 1

  with open('../ready_to_remodel/'+match,'w') as pdb_out:
    pdb_out.write(straight_As)
os.chdir('..')

### RENUMBER
os.chdir('ready_to_remodel/')
for remodelready in os.listdir('.'):
  if remodelready.endswith('.pdbA.pdb'):
    continue
  if remodelready.endswith('.pdb'):
    print 'will renumber: '+remodelready
    os.system('../clean_pdb.py '+remodelready+' A')
os.chdir('..')

### BLUEPRINTS
print "pose numbers:",temp_pose_num_dict
os.chdir('ready_to_remodel/')
for remodelready in os.listdir('.'):
  if remodelready.endswith('.pdbA.pdb'):
    this_id = remodelready[:4]
    blu_name = this_id + '.blueprint.template'
    os.system('../getBluePrintFromCoords.pl -pdbfile '+remodelready+' > '+blu_name)
    
    for this_length in range(3,9):
      new_blueprint = []
      with open(blu_name,'r') as blueprint_in:
        for line in blueprint_in.readlines():
          linelist = list(line)
          if int(line.split()[0]) == temp_pose_num_dict[this_id]['botA']:
            linelist[-2] = 'L'
            for x in range(this_length):
              linelist.append('0   x L ALLAA\n')
          if int(line.split()[0]) == temp_pose_num_dict[this_id]['capA']:
            linelist[-2] = 'L'
          if int(line.split()[0]) == temp_pose_num_dict[this_id]['capB']:
            linelist[-2] = 'L'
            for x in range(this_length):
              linelist.append('0   x L ALLAA\n')
          if int(line.split()[0]) == temp_pose_num_dict[this_id]['botB']:
            linelist[-2] = 'L'
          new_blueprint.append(''.join(linelist))
      with open( this_id+'.'+str(this_length)+'.blueprint', 'w') as blueprint_out:
        blueprint_out.write(''.join(new_blueprint))
os.chdir('..')

### REMODEL
print 'combos to remodel: ', temp_pose_num_dict.keys() #(hacky)
os.chdir('remodeled.fa.ALLAA')
bluelist = [ blue for blue in os.listdir('../ready_to_remodel/') ]
bluelist.sort()
for blueprint in bluelist:
  if blueprint.endswith('.blueprint'):
    this_id = blueprint[:4]
    combo_of_interest = this_id+'.trim.pdbA.pdb'
    print blueprint
    os.system( '../wrapper.sh '+combo_of_interest+' '+blueprint )
    
    if '.8.' in blueprint:
      print "this round complete" #raw_input('first combo done')

