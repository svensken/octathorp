#!/usr/bin/env python

from tempfile import mkstemp
from shutil import move
import os
from os import remove, close



params_template = """
hi

"""



def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)
    for line in old_file:
        new_file.write(line.replace(pattern, subst))
        #close temp file
        new_file.close()
        close(fh)
        old_file.close()
        #Remove original file
        remove(file_path)
        #Move new file
        move(abs_path, file_path)

def distance(x1, y1, z1, x2, y2, z2):
    dist = sqrt( (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )




bottoms = ['a']#['1WM1bot.pdb', '2WUFbot.pdb', '3ANSbot.pdb', '3B12bot.pdb', '3GZJbot.pdb']
caps = ['b']#['1WM1cap.pdb', '2WUFcap.pdb', '3ANScap.pdb', '3B12cap.pdb', '3GZJcap.pdb']

for bot in bottoms:
    for cap in caps:
        os.mkdir( bot+'_'+cap )
        os.chdir( bot+'_'+cap )
        with open('params.txt','w') as params:
            params.write( params_template )

        os.chdir( '..' )
