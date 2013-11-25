#!/usr/bin/env python

import os

dirr = 'scores/'

for scorefile in os.listdir(dirr):
  with open(dirr+scorefile, 'r') as inn:
    alllines=inn.readlines()
    t=[]
    for a in alllines:
      t.append( ','.join( a.split() ) )
  with open(dirr+scorefile+'.csv', 'w') as outt:
    outt.write('\n'.join(t))


