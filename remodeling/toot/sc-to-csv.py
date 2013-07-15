
with open('default.ssv', 'r') as inn:                                                                                               
  alllines=inn.readlines()                                                                                                          
  t=[]
  for a in alllines:                                                                                                                
    t.append( ','.join( a.split() ) )
    print a
  with open('out.csv', 'w') as outt:                                                                                                
    outt.write('\n'.join(t))                                                                                                        


