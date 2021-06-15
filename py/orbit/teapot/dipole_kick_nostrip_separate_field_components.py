import sys
import os
import math

# import teapot base functions from wrapper around C++ functions
from orbit.teapot_base import TPB

# import the function that creates multidimensional arrays
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import the MAD parser to construct lattices of TEAPOT elements.
from orbit.parsers.mad_parser import MAD_Parser, MAD_LattElement
# import the MADX parser to construct lattices of TEAPOT elements.
from orbit.parsers.madx_parser import MADX_Parser, MADX_LattElement

# import aperture
from aperture import Aperture

# monitor
from bunch import BunchTwissAnalysis

from orbit.teapot.teapot import BaseTEAPOT

from orbit_utils import Function

#computes the tracking functions and stripping probabilities
from orbit.teapot.function_strippingIncludeChicaneField import probabilityStrippingWithChicane

debugPrint=True
def printDebug(m1,m2="",m3="",m4="",m5="",m6="",m7="",m8="") :
	if debugPrint==True:
		print m1,m2,m3,m4,m5,m6,m7,m8
		
class GeneralDipoleNoStripSeperateField(BaseTEAPOT):
	def __init__(self,magneticFieldx,magneticFieldy,nParts,length,gamma,beta, name = "GeneralDipoleNoStripSeperateField no name"):
		"""
		Constructor. Creates the Dipole Combined Functions TEAPOT element .
		"""
		BaseTEAPOT.__init__(self, name)
		self.setType("parent stripper dipole no strip teapot")
		
		self.magneticFieldx=magneticFieldx
		self.magneticFieldy=magneticFieldy
		self.nParts=nParts
		self.length=length
		self.setLength(self.length)
		self.gamma=gamma
		self.beta=beta
			
		self.functionXPRigidity=None
		self.functionXRigidity=None
		
		self.functionYPRigidity=None
		self.functionYRigidity=None
		
		self.strippingTracking=None
		
		self.chicaneFieldx=0.
		self.chicaneFieldy=0.	
		
		self.computeFunctions()		


	def setChicaneFieldx(self,value):
		self.chicaneFieldx=value
	def getChicaneFieldx(self):
		return self.chicaneFieldx
	def setChicaneFieldy(self,value):
		self.chicaneFieldy=value
	def getChicaneFieldy(self):
		return self.chicaneFieldy		
	def setFunctionMagneticFieldx(self,function):
		self.magneticFieldx=function
	def getFunctionMagneticFieldx(self):
		return self.magneticFieldx			
	def setFunctionMagneticFieldy(self,function):
		self.magneticFieldy=function
	def getFunctionMagneticFieldy(self):
		return self.magneticFieldy	
		
	def setnParts(self,nParts):
		self.nParts=nParts
	def getnParts(self):
		return self.nParts	
	def setgamma(self,gamma):
		self.gamma=gamma
	def getgamma(self):
		return self.gamma
	def setbeta(self,beta):
		self.beta=beta
	def getbeta(self):	
		return self.beta	
		
	def setEffLength(self,length=0.01):
		#sets the effective length of the magnet
		self.length=length
		
	def getEffLength(self):
		#gets the effective length of the magnet
		return self.length
		
	def setFunctionCDF(self,function):
		self.functionCDF=function
	def getFunctionCDF(self):
		return self.functionCDF		

	def setFunctionXPRigidity(self,function):
		self.functionXPRigidity=function
	def getFunctionXPRigidity(self):
		return self.functionXPRigidity		
		
	def setFunctionXRigidity(self,function):
		self.functionXRigidity=function
	def getFunctionXRigidity(self):
		return self.functionXRigidity	
		
	def setFunctionYPRigidity(self,function):
		self.functionYPRigidity=function
	def getFunctionYPRigidity(self):
		return self.functionYPRigidity		
		
	def setFunctionYRigidity(self,function):
		self.functionYRigidity=function
	def getFunctionYRigidity(self):
		return self.functionYRigidity	
	def computeFunctions(self):
		self.strippingTracking=probabilityStrippingWithChicane(self.magneticFieldx,self.magneticFieldy,self.nParts,self.length,self.gamma,self.beta,False)
		self.strippingTracking.computeFunctions()		
		
		self.functionXPRigidity=self.strippingTracking.getdeltaxp_rigidity()
		self.functionXRigidity=self.strippingTracking.getdeltax_rigidity()

		self.functionYPRigidity=self.strippingTracking.getdeltayp_rigidity()
		self.functionYRigidity=self.strippingTracking.getdeltay_rigidity()

		self.functionInverse=self.strippingTracking.getInverseFunction()
				
	def track(self, paramsDict):
		"""
		The Dipole Combined Functions TEAPOT  class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getEffLength()

		bunch = paramsDict["bunch"]

		TPB.dipoleGeneralNoKickStripSeperateField(bunch,self.functionXPRigidity,self.functionXRigidity,self.functionYPRigidity,self.functionYRigidity,length)


	