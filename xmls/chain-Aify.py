#!/usr/bin/env python

import sys

filename = sys.argv[1]

with open(filename, 'r') as inn: 
  alllines=inn.readlines() 
  t=[]
  for a in alllines:
    if a.startswith('ATOM') and a[21] == 'B':
      a=list(a)
      a[21] = 'A'
      a=''.join(a)
      t.append( a )
    else:
      t.append( a )
    print a
  with open(filename, 'w') as outt:
    outt.write('\n'.join(t)) 


