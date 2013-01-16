import sys


try:
    pdb_file = [arg for arg in sys.argv if arg[-4:] == ".pdb"][0]
    print pdb_file

except:
    print "\nUsage: \n\
    $ draft.py [file.pdb] [int residue 1] [int residue 2] \n\
    \n\
    Please put all relevant PDB files in current directory.\n"
#    Please give PDB file path relative to current directory.


