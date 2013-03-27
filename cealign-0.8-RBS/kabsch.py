#!python

############################################################################
#
#  Copyright (c) 2007, Jason Vertrees.
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions
#  are
#  met:
#
#      * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#      * Redistributions in binary form must reproduce the above
#      copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the
#      distribution.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS
#  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
#  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
#  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  OWNER
#  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#  SPECIAL,
#  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
#  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
#  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
#  OF
#  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
#  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#############################################################################

##############################################################################
#
# @SUMMARY: -- Kabsch.py.  A python implementation of the optimal superposition
#     of two sets of vectors as proposed by Kabsch 1976 & 1978.
#
# @AUTHOR: Jason Vertrees
# @COPYRIGHT: Jason Vertrees (c), 2005.
# DATE  : 2005-04-07
# REV   : 2
# NOTES: Updated RMSD, notes, cleaned up the code a little.
#
#############################################################################
# math imports
import math
import Numeric
import LinearAlgebra
import Matrix

from array import *

# system stuff
import os
import copy

# pretty printing
import pprint
import string

# for importing as a plugin into PyMol
#import tkSimpleDialog
#import tkMessageBox
from pymol import cmd
from pymol import stored
from pymol import selector

class kabsch:
	"""
	Kabsch alignment of two set of vectors to produce and optimal alignemnt.

	Steps
	=====	
		1.  Calculate the center of mass for each protein.  Then, move the protein
		to its center of mass.  We choose as a convention, to use the origin as 
		the center of mass of the first set of coordinates.  This will allow us to
		return one translation vector, instead of two.
		
		Update: Since any rotation around a point that's not the origin, is in fact
		an affine rotation.  So, to beat that, we translate both to the center.
		
		NAME: superpose(c1, c2)
		POST: T is now defined.
		
		2.  Calculate the matrix, R, by (eq7, 1976).  r_{i,j} = sum_n w_n * y_{ni} * x_{nj},
		where y_{ni} is the ith component of the vector y_n.
		NAME: calcR
		POST: R is now defined
		
		3.  Calculate RtR (R-transpose * R).
		NAME: calcRtR
		POST: RtR is now defined
		
		4.  Calculate the corresponding eigenpairs for RtR, and sort them accordingly:
		m1 >= m2 >= m3; set v3 = v1 x v2 to ensure a RHS
		NAME: calcEigenPairs
		POST: The eigen pairs are calculated, sorted such that m1 >= m2 >= m3 and
		v3 = v1 x v2.
		
		5.  Calculate R*v_k and normalize the first two vectors to obtain b_1, b_2 and set
		b_3 = b_1 x b_2.  This also takes care of m1 > m2 = 0. 
		NAME: calcBVectors
		POST: b-Vectors are defined and returned
		
		6.  Calculate U=(u_{ij})=(sum n b_{ki} * a_{kj}) to obtain the best rotation.  Set
		sigma_3 = -1 if b3(Ra3) < 0 else sigma_3 = +1.
		NAME: calcU
		POST: U is defined
		
		7.  Calculate the RMSD.  The residual error is then
		The E = E0 - sqrt(m1) - sqrt(m2) - sigma_3(sqrt(m3)).
		NAME: calcRMSD
		POST: RMSD is computed.
	
	 @note: This should be a static method that takes three parameters
		
		1. The first protein's coordinates, f.  This program will 
		accept coordinates in the form an array of 3D vectors/lists
		
		2. The second protein's coordinates, g.
		
		3. The array of integer pairs representing the pairs to align.
		Coordinates should be formatted as as array of 2D vectors/lists.
	

	
	"""
	
	def __init__(self):
		"""
		Constructor.  @see kabsch.align.
		
		"""

		#
		# Instance Variables:  All three of these will be updated
		# every time the user calls ~.align.  Just to warn ya'.
		#		
		# U, the rotation matrix
		self.U = []
		# T, the translation vector
		self.T = []
		# R, the RMSD
		self.R = -1.0

		#self.menuBar.addmenuitem('Plugin', 'command', 'Kabsch Align', label = "Align Selections to Optial RMSD", command = lamda s=self: fetchPDBDialog(s))
		
				
	def align(self, c1, c2, pairs):
		"""
		Finds the best alignment of c1 and c2's pairs resulting in
		the smallest possible RMSD.

				
		@note:
			- All weights in this first version are set to 1.  Kabsch allows,
			differential weighting.  In the future, I may extend to this option,
			and then pairs may become 3D (r1, r2, weight) or I may add another
			parameter.
		
			- Helper functions will soon be provided such that the user may
			just use this package alone to compute the rotation.
		
		@param c1: coordinats of the first vectors, as an array of 3D vectors.
		@type  c1: Python list
		   
		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
		@type  c2: Python list
		   
		@param pairs: the list of pairs as an array of 2D pairs.
		@type  pairs: Python list of 2D lists.
		
		@return: U, the rotation matrix that gives the optimal rotation between the two proteins.
			
			T, the translation vector given to align the two objects centers of mass.
			
			R, the RMSD of the final rotation.
		"""
		
		#
		# First we move the center of mass of one protein, to the
		# center of mass of the other.  This removes any translation
		# between the two.
		#
		T1, T2, c1, c2 = self.superpose(c1, c2)
		# Calculate the initial RMSD
		E0 = self.calcE0(c1, c2)
		# Calculate R via eq. 7.
		R = self.calcR(c1, c2)
		# Calculate R(transpose)*R
		RtR = self.calcRtR(R)
		# Determined the eigenpairs for the matrix RtR.
		eValues, eVectors = self.calcEigenPairs(RtR)
		# Determine the bVectors as required
		bVectors = self.calcBVectors(R, eVectors)
		# Calculate the roation matrix
		U = self.calcU(eVectors, bVectors)
		# Calculate the final RMSD using U.
		RMSD = self.calcRMSD(E0, eValues, eVectors, bVectors, R, len(c1))
		
		return U, T1, T2, RMSD, c1, c2
		
		
	def superpose(self, c1, c2 ):
		"""
		Calculate the center of mass for each protein.  Then, move the protein
		to its center of mass.  We choose as a convention, to use the origin as 
		the center of mass of the first set of coordinates.  This will allow us to
		return one translation vector, instead of two.
		(CORRECT)
		
		@precondition: c1 and c2 are well defined lists of N-dimensional points with length > 0.
		@postcondition: T is now defined.
		
		@param c1: coordinats of the first vectors, as an array of 3D vectors.
		@type  c1: Python list
		   
		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
		@type  c2: Python list
		
		@return: T the translation vector.
		
			c2 one list of coordinates that of c2 translated to the COM of c1.
		
		"""
		
		# make sure we don't get bad data
		if (len(c1) != len(c2)):
			print "Two different length selections, with lengths, %d and %d." % (len(c1), len(c2))
			print "This algorithm must be used with selections of the same length."	
			print "In PyMol, type 'count_atoms sel1' where sel1 are your selections to find out their lengths."
			print "Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122"
        		print "Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA"

			
		assert(len(c1) == len(c2) != 0)

		L = len(c1)
		
		#
		# Centers of Mass
		#
		c1COM = Numeric.zeros((3,1), Numeric.Float64)
		c2COM = Numeric.zeros((3,1), Numeric.Float64)

		# calculate the CsOM		
		for i in range(0, L):
			for j in range(0,3):
				c1COM[j] = c1COM[j] + c1[i][j]
				c2COM[j] = c2COM[j] + c2[i][j]
		
		T1 = - c1COM / L
		T2 = - c2COM / L

		# move everything back to the origin.
		for i in range(0, L):
			for j in range(0,3):
				c1[i][j] = c1[i][j] + T1[j]
				c2[i][j] = c2[i][j] + T2[j]
				
		return T1, T2, c1, c2

					
	def calcR( self, c1, c2 ):
		"""
		Calculate the matrix, R, by (eq7, 1976).  M{r_{i,j} = sum_n w_n * y_{ni} * x_{nj}},
		where M{y_{ni}} is the ith component of the vector M{y_n}.
		(CORRECT)
		
		@param c1: coordinats of the first vectors, as an array of 3D vectors.
		@type  c1: Python list
		   
		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
		@type  c2: Python list

		@postcondition: R is now defined.
		
		@return: R
		@rtype : 3x3 matrix
				
		"""
		
		# Create the 3x3 matrix
		R = Numeric.zeros((3,3), Numeric.Float64)
		L = len(c1)
		
		for k in range(0, L):
			for i in range(0, 3):
				for j in range(0, 3):
					R[i][j] = R[i][j] + (c2[k][i] * c1[k][j])
		
		# return R the 3x3 PSD Matrix.
		return R
		
	
	def calcRtR( self, R ):
		"""
		Calculate RtR (R-transpose * R).
		(CORRECT)
		
		@param R: Matrix
		@type  R: 3x3 Matrix
		
		@precondition: R is a the well formed matrix as per Kabsch.
		@postcondition: RtR is now defined
		
		@return: M{R^tR}
		@rtype : 3x3 matrix
		
		"""
		
		RtR = Numeric.matrixmultiply(Numeric.transpose(R), R)
		
		return RtR
	
	
	def calcEigenPairs( self, RtR ):
		"""
		Calculate the corresponding eigenpairs for RtR, and sort them accordingly:
		M{m1 >= m2 >= m3}; set M{v3 = v1 x v2} to ensure a RHS
		(CORRECT)
		
		@postcondition: The eigen pairs are calculated, sorted such that M{m1 >= m2 >= m3} and
		M{v3 = v1 x v2}.
		
		@param RtR: 3x3 Matrix of M{R^t * R}.
		@type  RtR: 3x3 Matrix
		@return : Eigenpairs for the RtR matrix.
		@rtype  : List of stuff
		
		"""
		
		eVal, eVec = LinearAlgebra.eigenvectors(RtR)

		# This is cool.  We sort it using Numeric.sort(eVal)
		# then we reverse it using nifty-crazy ass notation [::-1].
		eVal2 = Numeric.sort(eVal)[::-1]
		eVec2 = [[],[],[]] #Numeric.zeros((3,3), Numeric.Float64)
				
		# Map the vectors to their appropriate owners		
		if ( eVal2[0] == eVal[0]):
			eVec2[0] = eVec[0]
			if ( eVal2[1] == eVal[1] ):
				eVec2[1] = eVec[1]
				eVec2[2] = eVec[2]
			else:
				eVec2[1] = eVec[2]
				eVec2[2] = eVec[1]
		elif( eVal2[0] == eVal[1]):
			eVec2[0] = eVec[1]
			if ( eVal2[1] == eVal[0] ):
				eVec2[1] = eVec[0]
				eVec2[2] = eVec[2]
			else:
				eVec2[1] = eVec[2]
				eVec2[2] = eVec[0]
		elif( eVal2[0] == eVal[2]):
			eVec2[0] = eVec[2]
			if ( eVal2[1] == eVal[1] ):
				eVec2[1] = eVec[1]
				eVec2[2] = eVec[0]
			else:
				eVec2[1] = eVec[0]
				eVec2[2] = eVec[1]

		eVec2[2][0] = eVec2[0][1]*eVec2[1][2] - eVec2[0][2]*eVec2[1][1]
		eVec2[2][1] = eVec2[0][2]*eVec2[1][0] - eVec2[0][0]*eVec2[1][2]
		eVec2[2][2] = eVec2[0][0]*eVec2[1][1] - eVec2[0][1]*eVec2[1][0]
		
		return [eVal2, eVec2]
		
	
	def calcBVectors( self, R, eVectors ):
		"""
		Calculate M{R*a_k} and normalize the first two vectors to obtain M{b_1, b_2} and set
		M{b_3 = b_1 x b_2}.  This also takes care of {m2 > m3 = 0}. 
		(CORRECT)
		
		@postcondition: b-Vectors are defined and returned
		
		@return: The three B-vectors
		@rtype: List of 3D vectors (Python LOL).
		"""
		
		bVectors = Numeric.zeros((3,3), Numeric.Float64)

		for i in range(0,3):
			bVectors[i] = Numeric.matrixmultiply(R, eVectors[i])

		bVectors[0] = bVectors[0] / Numeric.sqrt(Numeric.add.reduce(bVectors[0]**2))
		bVectors[1] = bVectors[1] / Numeric.sqrt(Numeric.add.reduce(bVectors[1]**2))
		bVectors[2] = bVectors[2] / Numeric.sqrt(Numeric.add.reduce(bVectors[2]**2))
		
		bVectors[2][0] = bVectors[0][1]*bVectors[1][2] - bVectors[0][2]*bVectors[1][1]
		bVectors[2][1] = bVectors[0][2]*bVectors[1][0] - bVectors[0][0]*bVectors[1][2]
		bVectors[2][2] = bVectors[0][0]*bVectors[1][1] - bVectors[0][1]*bVectors[1][0]
		
		return bVectors
		
		
		
	def calcU( self, eVectors, bVectors ):
		"""
		Calculate M{U=(u_{ij})=(sum n b_{ki} * a_{kj})} to obtain the best rotation.  Set
		M{sigma_3 = -1 if b3(Ra3) < 0 else sigma_3 = +1}.
		(CORRECT)
		
		@postcondition: U is defined
		
		@param eVectors: Eigenvectors for the system.
		@type  eVectors: Eigenvectors
		
		@param bVectors: BVectors as described by Kabsch.
		@type  bVectors: BVectors
		
		@return: U the rotation matrix.
		@rtype  :3x3 matrix.
		
		"""
		
		U = Numeric.zeros((3,3), Numeric.Float64)
		
		for k in range(0,3):
			for i in range(0,3):
				for j in range(0,3):
					U[i][j] = U[i][j] + Numeric.matrixmultiply(bVectors[k][i], eVectors[k][j])
		
		return U
	
		
	def calcE0( self, c1, c2 ):
		"""
		Calculates the initial RMSD, which Kacbsch called E0.
		(CORRECT)
		
		@param c1: coordinats of the first vectors, as an array of 3D vectors.
		@type  c1: Python list
		   
		@param c2: coordinates of the second set of vectors, as an array of 3D vectors.
		@type  c2: Python list
		
		@return: E0 the initial RMSD.
		@rtype : float.
				
		"""
		
		E0 = 0.0
		
		L = len(c1)
		for i in range( 0, L ):
			for j in range(0, 3):
				E0 = E0 + 0.5*( (c1[i][j]*c1[i][j])+(c2[i][j]*c2[i][j]))
		
		return E0
	
	def calcRMSD( self, E0, eValues, eVectors, bVectors, R, N):
		"""
		Calculate the RMSD.  The residual error is then
		The M{E = E0 - sqrt(m1) - sqrt(m2) - sigma_3(sqrt(m3))}.
		
		@param E0: Initial RMSD as calculated in L{calcE0}.
		@type  E0: float.
		
		@param eVectors: Eigenvectors as calculated from L{calcEigenPairs}
		@type eVectors: vectors, dammit!
		
		@param bVectors: B-vectors as calc. from L{calcBVectors}
		@type  bVectors: More vectors.
		
		@param R: The matrix R, from L{calcR}.
		@type  R: 3x3 matrix.		

		@param N: Number of equivalenced points
		@type  N: integer
		
		@postcondition: RMSD is computed.
		@return: The RMSD.
		
		"""
		sigma3 = 0
		if ( Numeric.matrixmultiply(bVectors[2], Numeric.matrixmultiply( R, eVectors[2])) < 0):
			sigma3 = -1
		else:
			sigma3 = 1

		E = Numeric.sqrt( 2*(E0 - Numeric.sqrt(eValues[0]) - Numeric.sqrt(eValues[1]) - sigma3*(Numeric.sqrt(eValues[2]))) / N)
		
		return E
		
		
	def calcSimpleRMSD( self, c1, c2 ):
		"""
		Calculates the usual concept of RMSD between two set of points.  The CalcRMSD above
		sticks to Kabsch's alignment method protocol and calculates something much different.
		@see kabsch.calcRMSD
		
		@param c1: List of points #1
		@type  c1: LOL
		
		@param c2: List of points #2
		@type  c2: LOL
		
		@return: RMSD between the two
		
		"""
		
		RMSD = 0.0
		for i in range(0, len(c1)):
			for j in range(0,3):
				RMSD = RMSD + (c2[i][j]-c1[i][j])**2
				
		RMSD = RMSD / len(c1)
		RMSD = Numeric.sqrt(RMSD)
		return RMSD
		
		
	#####################################################################
	#
	# UTIL Functions
	def rotatePoints(self, U, c2):
		"""
		Rotate all points in c2 based on the rotation matrix U.
		
		@param U: 3x3 Rotation matrix
		@type  U: 3x3 matrix
		
		@param c2: List of points to rotate
		@type  c2: List of 3D vectors
		
		@return: List of rotated points
		
		"""
	
		L = len(c2)

		for n in range(0,L):
			c2[n][0] = c2[n][0] * U[0][0] + c2[n][1] * U[1][0] + c2[n][2] * U[2][0]
			c2[n][1] = c2[n][0] * U[0][1] + c2[n][1] * U[1][1] + c2[n][2] * U[2][1]
			c2[n][2] = c2[n][0] * U[0][2] + c2[n][1] * U[1][2] + c2[n][2] * U[2][2]
		
		return c2
	
	def writeU( self, U, fileName ):
		"""
		Convenience function.  Writes U to disk.
		
		"""
		
		if ( len(fileName) == 0 ):
			fileName = "./U"
			
		outFile = open( fileName, "wb")
		for i in range(0,3):
			for j in range(0,3):
				outFile.write( string.ljust(str(U[i][j]),20) )
			outFile.write("\n")
		outFile.close()		

	

def optAlign( sel1, sel2 ):
	"""
	optAlign performs the Kabsch alignment algorithm upon the alpha-carbons of two selections.
	Example:   optAlign MOL1 and i. 20-40, MOL2 and i. 102-122
	Example 2: optAlign 1GGZ and i. 4-146 and n. CA, 1CLL and i. 4-146 and n. CA
	
	Two RMSDs are returned.  One comes from the Kabsch algorithm and the other from
	PyMol based upon your selections.

	By default, this program will optimally align the ALPHA CARBONS of the selections provided.
	To turn off this feature remove the lines between the commented "REMOVE ALPHA CARBONS" below.
	
	@param sel1: First PyMol selection with N-atoms
	@param sel2: Second PyMol selection with N-atoms
	"""
	cmd.reset()

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
	sel1 = sel1 + " and N. CA"
	sel2 = sel2 + " and N. CA"
	# -- REMOVE ALPHA CARBONS

	cmd.iterate_state(1, selector.process(sel1), "stored.sel1.append([x,y,z])")
	cmd.iterate_state(1, selector.process(sel2), "stored.sel2.append([x,y,z])")
	# full molecule
	mol1 = cmd.identify(sel1,1)[0][0]
	mol2 = cmd.identify(sel2,1)[0][0]
	cmd.iterate_state(1, mol1, "stored.mol1.append([x,y,z])")
	cmd.iterate_state(1, mol2, "stored.mol2.append([x,y,z])")

	K = kabsch()
	U, T1, T2, RMSD, c1, c2 = K.align(stored.sel1, stored.sel2, [])

	stored.mol2 = map(lambda v:[T2[0]+((v[0]*U[0][0])+(v[1]*U[1][0])+(v[2]*U[2][0])),T2[1]+((v[0]*U[0][1])+(v[1]*U[1][1])+(v[2]*U[2][1])),T2[2]+((v[0]*U[0][2])+(v[1]*U[1][2])+(v[2]*U[2][2]))],stored.mol2)
	#stored.mol1 = map(lambda v:[ v[0]+T1[0], v[1]+T1[1], v[2]+T1[2] ], stored.mol1)
	stored.mol1 = map(lambda v:[ v[0]+T1[0], v[1]+T1[1], v[2]+T1[2] ], stored.mol1)

	cmd.alter_state(1,mol1,"(x,y,z)=stored.mol1.pop(0)")
	cmd.alter_state(1,mol2,"(x,y,z)=stored.mol2.pop(0)")
	cmd.alter( 'all',"segi=''")
	cmd.alter('all', "chain=''")
	print "RMSD=%f" % cmd.rms_cur(sel1, sel2)
	print "MY RMSD=%f" % RMSD
	cmd.hide('everything')
	cmd.show('ribbon', sel1 + ' or ' + sel2)
	cmd.color('gray70', mol1 )
	cmd.color('paleyellow', mol2 )
	cmd.color('red', 'visible')
	cmd.show('ribbon', 'not visible')
	cmd.center('visible')
	cmd.orient()
	cmd.zoom('visible')

cmd.extend("optAlign", optAlign)	

