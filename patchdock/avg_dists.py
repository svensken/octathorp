#!/usr/bin/env python

import os, math



def distance(x1, y1, z1, x2, y2, z2):
    dist = math.sqrt( (float(x1)-float(x2))**2 + (float(y1)-float(y2))**2 + (float(z1)-float(z2))**2 )
    return dist



p = '/home/svensken/octathorp/patchdock/'

#bottoms = {'1WM1bot.pdb':'138_248', '2WUFbot.pdb':'139_222', '3ANSbot.pdb':'361_485', '3B12bot.pdb':'128_230', '3GZJbot.pdb':'111_198'}
#caps = {'1WM1cap.pdb':'140_238', '2WUFcap.pdb':'145_216', '3ANScap.pdb':'369_474', '3B12cap.pdb':'139_226', '3GZJcap.pdb':'121_195'}

# for lu
bottoms = {'1WM1':'138_248', '2WUF':'139_222', '3ANS':'361_485', '3B12':'128_230', '3GZJ':'111_198'}
caps = {'1WM1':'140_238', '2WUF':'145_216', '3ANS':'369_474', '3B12':'139_226', '3GZJ':'121_195'}


#for bot in bottoms:
#    for cap in caps:

# for lu
for bot in [ bottoms['1WM1'] ]:
    for cap in [ caps['1WM1'] ]:

        if bot[:4] == cap[:4]:
            # assuming cealigned
            print "NATIVE :D"

        #os.chdir( 'combos/'+bot+'_'+cap )
        
        # for lu
        os.chdir( 'combos/'+bot+'_'+cap )


        with open('dists.csv','w') as d:
            d.write('dir,dist1,dist2\n') 

        for pdb in os.listdir('.'):
            if pdb.startswith('output') and pdb.endswith('.pdb'):

                with open( pdb, 'r' ) as pdbfile:
                    lines = pdbfile.readlines()

                    bot_res_nums = bottoms[bot].split('_')
                    bot_res1 = bot_res_nums[0]
                    bot_res2 = bot_res_nums[1]
                    cap_res_nums = caps[cap].split('_')
                    cap_res1 = cap_res_nums[0]
                    cap_res2 = cap_res_nums[1]

                    for line in lines:
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
                            pass

                dist1 = distance( bot_x1, bot_y1, bot_z1,
                                  cap_x1, cap_y1, cap_z1 )
                dist2 = distance( bot_x2, bot_y2, bot_z2,
                                  cap_x2, cap_y2, cap_z2 )

                with open('dists.csv','a') as d:
                    linetowrite = bot+'_'+cap+','+str(dist1)+','+str(dist2)
                    d.write( linetowrite+'\n' )

        os.chdir( p )
