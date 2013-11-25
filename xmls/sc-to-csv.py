#!/usr/bin/env python

import os

scorefile = '3ANS_3ANS.sc'

with open(scorefile, 'r') as inn:
  alllines=inn.readlines()
  t=[]
  for a in alllines:
    t.append( ','.join( a.split() ) )
with open(scorefile+'.csv', 'w') as outt:
  outt.write('\n'.join(t))


