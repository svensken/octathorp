#
# -- cealign.py
#

############################################################################
#
#  Copyright (c) 2007, Jason Vertrees.
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are
#  met:
#  
#      * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#      * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the
#      distribution.
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
#  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
#  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
#  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#############################################################################

#############################################################################
#
# Purpose: Provides a PyMol plugin that calculates the 3D structural align-
# ======== ment for two proteins.
#
# Notes  : This performs the first step of the CE alignment algorithm.  It
# ======== does not perform the final step of optimization of alignment for
#          high-scoring sequences.
#
# Installation: See provided documentation (READMEs, INSTALLs, and Python
# ============= docstrings.  Briefly, you'll need: Kabsch, or QKabsch,
#               from the wiki.  Also, you'll need to install numpy.
#
# Author : Jason Vertrees.
# ========
#
#############################################################################
from pymol import cmd, selector, stored
import pprint
import math
import copy


from ccealign import ccealign

def cealign( sel1, sel2 ):
	"""
	Rough CE Structure-based Alignment of two protein structures

	Overview:
	cealign will try its best to align the alpha-carbon atoms provided
	in both selections sel1 and sel2.  The algorithm follows the paper
	written by Shindyalov and Bourne.  A few modifications to the algo-
	rithm are introduced, partly due to lazyness, and partly due to 
	improving speed and accuracy while not conflicting with the lazyness
	requirement.  :-)

	Params:
	\@param sel1: (string) valid PyMol selection string of protein 1 to align	
	\@param sel2: (string) valid PyMol selection string of protein 2 to align

	Returns:
	\@return: (string) the CE-score and the RMSD of the alignment

	Side-Effects:
	\@note: rotates and translates the objects (proteins) provided in the
	selections, sel1 and sel2, to represent the alignment.  Probably will
	also change their representation to more clearly show the aligned
	segments.

	Requires:
	Requires the Kabsch algorithm for protein optimal superposition.  Don't
	worry, I already wrote this as a foray into the academic: you may download
	and install it from the PyMol wiki at:
		http://www.pymolwiki.org/index.php/Kabsch

	Bugs: Many, I'm sure.  Please forward bugs/comments to javertre@utmb.edu
	"""
	
	#########################################################################
	# CE specific defines.  Don't change these unless you know
	# what you're doing.  See the documentation.
	#########################################################################
	# WINDOW SIZE
	# make sure you set this variable in cealign.py, as well!
	winSize = 8
	# FOR AVERAGING
	winSum = (winSize-1)*(winSize-2) / 2;
	# max gap size
	gapMax = 30
	
	
	# make the lists for holding coordinates
	# partial lists
	stored.sel1 = []
	stored.sel2 = []
	# full lists
	stored.mol1 = []
	stored.mol2 = []
 
	# now put the coordinates into a list
	# partials
 
	# -- REMOVE ALPHA CARBONS
	sel1 = sel1 + " and n. CA"
	sel2 = sel2 + " and n. CA"
	# -- REMOVE ALPHA CARBONS
 
	cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
	cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
	
	# full molecule
	mol1 = cmd.identify(sel1,1)[0][0]
	mol2 = cmd.identify(sel2,1)[0][0]
	
	# put all atoms from MOL1 & MOL2 into stored.mol1
	cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
	cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")
	
	if ( len(stored.mol1) == 0 ):
		print "ERROR: Your first selection was empty."
		return
	if ( len(stored.mol2) == 0 ):
		print "ERROR: Your second selection was empty."
		return
	
	# call the C function
	alignString = ccealign( (stored.sel1, stored.sel2) )
	
	if ( len(alignString) == 1 ):
		if ( len(alignString[0]) == 0 ):
			print "\n\nERROR: There was a problem with CEAlign's C Module.  The return value was blank."
			print "ERROR: This is obviously bad.  Please inform a CEAlign developer.\n\n"
			return
	
	bestPathID = -1
	bestPathScore = 100000
	bestStr1 = ""
	bestStr2 = ""
	
	# for each of the 20 possible alignments returned
	# we check each one for the best CE-Score and keep
	# that one.  The return val of ccealign is a list
	# of lists of pairs.
	for curAlignment in alignString:
		seqCount = len(curAlignment)
		matA = None
		matB = None
		
		if ( seqCount == 0 ):
			continue;
		
		for AFP in curAlignment:
			first, second = AFP
			if ( matA == None and matB == None ):
				matA = [ stored.sel1[first-1] ]
				matB = [ stored.sel2[second-1] ]
			else:
				matA.append( stored.sel1[first-1] )
				matB.append( stored.sel2[second-1] )
		
		curScore = simpAlign( matA, matB, mol1, mol2, stored.mol1, stored.mol2, align=0, L=len(matA) )
		
		#########################################################################
		# if you want the best RMSD, not CE Score uncomment here down
		#########################################################################
		#if ( curScore < bestPathScore ):
			#bestPathScore = curScore
			#bestMatA = matA
			#bestMatB = matB
		#########################################################################
		# if you want the best RMSD, not CE Score uncomment here up
		#########################################################################
		
		#########################################################################
		# if you want a proven, "better" alignment use the CE-score instead
		# Uncomment here down for CE-Score
		#########################################################################
		internalGaps = 0.0;
		for g in range(0, seqCount-1):
			if (not curAlignment[g][0] + 1 == curAlignment[g+1][0]):
				internalGaps += curAlignment[g+1][0]
			if ( not curAlignment[g][1] + 1 == curAlignment[g+1][1] ):
				internalGaps += curAlignment[g+1][1]
			
			aliLen = float( len(curAlignment))
			numGap = internalGaps;
			curScore = float((curScore/aliLen)*(1.0+(numGap/aliLen)));
		
		if ( curScore < bestPathScore ):
			bestPathScore = curScore
			bestMatA = matA
			bestMatB = matB
		#########################################################################
		# if you want a proven, "better" alignment use the CE-score instead
		# Uncomment here UP for CE-Score
		#########################################################################

	# align the best one string
	simpAlign(bestMatA, bestMatB, mol1, mol2, stored.mol1, stored.mol2, align=1, L=len(bestMatA))


def alignto(sel1):
	"""Just a quick & dirty multiple structure alignment"""
	for x in cmd.get_names("*"): 
		print "Aligning %s to %s" % (x, cmd.get_names(sel1)[0])
		cealign( sel1, x )
	


## Let PyMol know about the command
cmd.extend("cealign", cealign)

## Let PyMOL know about the alignto command
cmd.extend("alignto", alignto)
