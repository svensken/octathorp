#!/usr/bin/env python

import os, math

def distance(x1, y1, z1, x2, y2, z2):
  dist = math.sqrt( (float(x1)-float(x2))**2 + (float(y1)-float(y2))**2 + (float(z1)-float(z2))**2 )
  return dist

os.chdir('completed/')

for combo_dir in os.listdir('../combos/'):
  print "###"+combo_dir+"###"

  for combo_file in os.listdir('../combos/'+combo_dir):
    print "#"+combo_file+"#"
    with open('../combos/'+combo_dir+'/'+combo_file,'r') as file_in:
      file_in = file_in.readlines()

      x1 = 0
      for line in file_in:
        elements = line.split()

        if elements[0]=="ATOM":
          if elements[2]=="N":
            if elements[4] == "A" or elements[4] == "B":
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
                if len(str(elements[5])) > 3:
                  print elements
                else:
                  print "res_"+elements[5]
