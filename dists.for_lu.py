#!/usr/bin/env python

import os, sys, math



def distance(x1, y1, z1, x2, y2, z2):
    dist = math.sqrt( (float(x1)-float(x2))**2 + (float(y1)-float(y2))**2 + (float(z1)-float(z2))**2 )
    return dist



p = '/home/svensken/octathorp/patchdock/'

#bottoms = {'1WM1bot.pdb':'138_248', '2WUFbot.pdb':'139_222', '3ANSbot.pdb':'361_485', '3B12bot.pdb':'128_230', '3GZJbot.pdb':'111_198'}
#caps = {'1WM1cap.pdb':'140_238', '2WUFcap.pdb':'145_216', '3ANScap.pdb':'369_474', '3B12cap.pdb':'139_226', '3GZJcap.pdb':'121_195'}

# for lu
#bottoms = {'1WM1':'138_248', '2WUF':'139_222', '3ANS':'361_485', '3B12':'128_230', '3GZJ':'111_198'}
#caps = {'1WM1':'140_238', '2WUF':'145_216', '3ANS':'369_474', '3B12':'139_226', '3GZJ':'121_195'}


combo = '1WM1_3ANS'
bot_res1 = "140"
bot_res2 = "246"
cap_res1 = "370"
cap_res2 = "469"


#for bot in bottoms:
#    for cap in caps:

# for lu
#for bot in [ bottoms['1WM1'] ]:
#    for cap in [ caps['1WM1'] ]:


os.chdir(combo)


# filename is 1WM1bot_3ANScap.pdb.dists.csv
with open('dists.csv','w') as d:
    d.write('struct_name,dist1,dist2\n') 


for pdb in os.listdir('.'):
    print pdb
    
    if pdb.endswith('.pdb'):

        with open( pdb, 'r' ) as pdbfile:
            lines = pdbfile.readlines()

            for line in lines:
                if line.split()[0] != "ATOM":
                    continue
                try:
                    if line.split()[5]==bot_res1 and line.split()[2]=='CA':
                        bot_atomid1 = line.split()[1]
                        bot_x1 = line.split()[6]
                        bot_y1 = line.split()[7]
                        bot_z1 = line.split()[8]
                    if line.split()[5]==bot_res2 and line.split()[2]=='CA':
                        bot_atomid2 = line.split()[1]
                        bot_x2 = line.split()[6]
                        bot_y2 = line.split()[7]
                        bot_z2 = line.split()[8]
                    if line.split()[5]==cap_res1 and line.split()[2]=='CA':
                        cap_atomid1 = line.split()[1]
                        cap_x1 = line.split()[6]
                        cap_y1 = line.split()[7]
                        cap_z1 = line.split()[8]
                    if line.split()[5]==cap_res2 and line.split()[2]=='CA':
                        cap_atomid2 = line.split()[1]
                        cap_x2 = line.split()[6]
                        cap_y2 = line.split()[7]
                        cap_z2 = line.split()[8]
                except:
                    print line

        dist1 = distance( bot_x1, bot_y1, bot_z1,
                          cap_x1, cap_y1, cap_z1 )
        dist2 = distance( bot_x2, bot_y2, bot_z2,
                          cap_x2, cap_y2, cap_z2 )

        with open('dists.csv','a') as d:
            linetowrite = pdb+','+str(dist1)+','+str(dist2)
            d.write( linetowrite+'\n' )

