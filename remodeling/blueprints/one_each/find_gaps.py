#!/usr/bin/env python

import os, math

def distance(x1, y1, z1, x2, y2, z2):
  dist = math.sqrt( (float(x1)-float(x2))**2 + (float(y1)-float(y2))**2 + (float(z1)-float(z2))**2 )
  return dist


for f in os.listdir('.'):
  x1 = 0
  if f.endswith('.pdb'):
    print f
    with open(f, 'r') as pdb:
      lines = pdb.readlines()
      for line in lines:
        elements = line.split()

        if elements[0]=="ATOM":
          if elements[2]=="N":
            if x1:
              x1_prev = x1
              y1_prev = y1
              z1_prev = z1
            else:
              x1_prev = float(elements[6])
              y1_prev = float(elements[7])
              z1_prev = float(elements[8])
            x1 = float(elements[6])
            y1 = float(elements[7])
            z1 = float(elements[8])
            
            if distance(x1,y1,z1,x1_prev,y1_prev,z1_prev) > 4:
              print elements[5]
