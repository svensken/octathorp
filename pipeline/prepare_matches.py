#!/usr/bin/env python

import os, sys


temp_resnum_dict = { '1WM1botA' : 138,
                     '1WM1capA' : 140,
                     '1WM1capB' : 238,
                     '1WM1botB' : 248,
                     '2WUFbotA' : 139,
                     '2WUFcapA' : 145,
                     '2WUFcapB' : 216,
                     '2WUFbotB' : 222,
                     '3ANSbotA' : 361,
                     '3ANScapA' : 369,
                     '3ANScapB' : 474,
                     '3ANSbotB' : 485,
                     '3B12botA' : 128,
                     '3B12capA' : 139,
                     '3B12capB' : 226,
                     '3B12botB' : 230,
                     '3GZJbotA' : 111,
                     '3GZJcapA' : 121,
                     '3GZJcapB' : 195,
                     '3GZJbotB' : 198 }


os.chdir('matched')
for match in os.listdir('.'):
  print match
  this_id = match[:4]
  
  straight_As = []
  chunk_o_Bs = []
  with open(match, 'r') as pdb_in:
    lines = pdb_in.readlines()
    for line in lines:
      linelist = list(line)
      if line[:4] == "ATOM":
        if linelist[21]=='A' and line[12:15]==' N ' and int(line[23:26])==temp_resnum_dict[this_id+'botB']:
          #print 'prev: ',prev_resi
          #print 'resi: ',resi
          #print 'delta: ', delta
          straight_As.append('\nrepwace_me\n')
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
  
  with open('../ready_to_remodel/'+match,'w') as pdb_out:
    pdb_out.write(straight_As)
