#!/usr/bin/env python

import os
ls=[]
for pdb in os.listdir('.'):
    if not pdb.startswith('UM'):
        continue
    
    with open(pdb, 'r') as pdbl:
        lines = pdbl.readlines()
        for line in lines:
            liney = line.split()
            if liney[0] == 'pose':
                total = liney[-1]
            if liney[0] == 'rmsA':
                rmsA = liney[-1]
            if liney[0] == 'rmsB':
                rmsB = liney[-1]
            if liney[0] == 'rmstot':
                rmstot = liney[-1]

    ls.append( pdb +', '+ rmstot +', '+ rmsA +', '+ rmsB +', '+ total )

with open('my.csv','w') as csv:
    csv.write('\n'.join(ls))

print 'done'
