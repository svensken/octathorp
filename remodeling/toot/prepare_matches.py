#!/usr/bin/env python

import os, sys

for d in os.listdir('.'):
  if '.' not in d:
  #if d == sys.argv[1]:
    if d == 'outp' or d == 'outty':
      continue
      # all i did was `./.py >> outp` and now it keeps showing up
    print d
    os.chdir(d)
    
    for f in os.listdir('.'):
      #print f
      straight_As = []
      chunk_o_Bs = []
      # test
      #if 'UM_1_D335H524D496Y155Y157_output_3ANS_1.pdb' not in f:
      #  continue
      with open(f,'r') as inn:
        lines = inn.readlines()
        for line in lines:
          linelist = list(line)
          if line[:4] == "ATOM":
            try:
              prev_resi = resi
            except:
              prev_resi = int(line[23:26])
            resi = int(line[23:26])
            delta = resi - prev_resi
            

            if linelist[21] == 'A' and delta > 10:
              #if (prev_resi==138 and resi==284) or \
              #   (prev_resi==139 and resi==222) or \
              #   (prev_resi==361 and resi==485) or \
              #   (prev_resi==128 and resi==230) or \
              #   (prev_resi==111 and resi==198):
                   
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

      with open('/home/svensken/octathorp/remodeling/matched/'+d+'_'+f,'w') as outt:
        outt.write(straight_As)
    os.chdir('..')
