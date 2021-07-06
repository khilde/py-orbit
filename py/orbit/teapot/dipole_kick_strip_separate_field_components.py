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
		
class GeneralDipoleStripSeperateField(BaseTEAPOT):
	def __init__(self,magneticFieldx,magneticFieldy,nParts,length,gamma,beta, name = "GeneralDipoleStripSeperateField no name",stripLength=-1):
		"""
		Constructor. Creates the Dipole Combined Functions TEAPOT element .
		"""
		BaseTEAPOT.__init__(self,name)
		self.setType("parent stripper dipole teapot")
		
		self.magneticFieldx=magneticFieldx
		self.magneticFieldy=magneticFieldy
		self.nParts=nParts
		self.length=length
		self.setLength(self.length)
		self.gamma=gamma
		self.beta=beta
		self.stripLength=stripLength
		
		self.functionInverse=None
		self.functionCDF=None
		
		#tracking functions for a charged particle that is stripped to neutral
		self.functionXPRigidity=None
		self.functionXRigidity=None
		self.functionYPRigidity=None
		self.functionYRigidity=None		
		
		#tracking functions for a neutral particle that is stripped to have charge
		self.functionXP_mRigidity=None
		self.functionX_mRigidity=None		
		self.functionYP_mRigidity=None
		self.functionY_mRigidity=None	
		
		
		self.strippingTracking=None
		self.computeFunctions()
		
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

	def setFunctionInverse(self,function):
		self.functionInverse=function
	def getFunctionInverse(self):
		return self.functionInverse	
		
	def setFunctionXPRigidity(self,function):
		self.functionXPRigidity=function
	def getFunctionXPRigidity(self):
		return self.functionXPRigidity		
		
	def setFunctionXRigidity(self,function):
		self.functionXRigidity=function
	def getFunctionXRigidity(self):
		return self.functionXRigidity	
		
	def setFunctionXP_mRigidity(self,function):
		self.functionXP_mRigidity=function
	def getFunctionXP_mRigidity(self):
		return self.functionXP_mRigidity		
		
	def setFunctionX_mRigidity(self,function):
		self.functionX_mRigidity=function
	def getFunctionX_mRigidity(self):
		return self.functionX_mRigidity	
		
		
	def setFunctionYPRigidity(self,function):
		self.functionYPRigidity=function
	def getFunctionYPRigidity(self):
		return self.functionYPRigidity		
		
	def setFunctionYRigidity(self,function):
		self.functionYRigidity=function
	def getFunctionYRigidity(self):
		return self.functionYRigidity	
		
	def setFunctionYP_mRigidity(self,function):
		self.functionYP_mRigidity=function
	def getFunctionYP_mRigidity(self):
		return self.functionYP_mRigidity		
		
	def setFunctionY_mRigidity(self,function):
		self.functionY_mRigidity=function
	def getFunctionY_mRigidity(self):
		return self.functionY_mRigidity	
	def computeFunctions(self):
		self.strippingTracking=probabilityStrippingWithChicane(self.magneticFieldx,self.magneticFieldy,self.nParts,self.length,self.gamma,self.beta)
		self.strippingTracking.computeFunctions()
		
		self.functionCDF=self.strippingTracking.getCDF()
		
		self.functionXPRigidity=self.strippingTracking.getdeltaxp_rigidity()
		self.functionXRigidity=self.strippingTracking.getdeltax_rigidity()
		self.functionXP_mRigidity=self.strippingTracking.getdeltaxp_m_rigidity()
		self.functionX_mRigidity=self.strippingTracking.getdeltax_m_rigidity()
		
		self.functionYPRigidity=self.strippingTracking.getdeltayp_rigidity()
		self.functionYRigidity=self.strippingTracking.getdeltay_rigidity()
		self.functionYP_mRigidity=self.strippingTracking.getdeltayp_m_rigidity()
		self.functionY_mRigidity=self.strippingTracking.getdeltay_m_rigidity()	
		
		self.functionInverse=self.strippingTracking.getInverseFunction()
		
	def track(self, paramsDict):
		"""
		The Dipole Combined Functions TEAPOT  class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getEffLength()

		bunch = paramsDict["bunch"]
		firstChicaneFail=paramsDict["firstChicaneFail"]
		secondChicaneFail=paramsDict["secondChicaneFail"]
		
		#charge <0 means first dipole
		#charge==0 means second dipole
		#TPB.dipoleGeneralKickStripSeperateField(bunch,firstChicaneFail,self.functionCDF,self.functionInverse,self.functionXPRigidity,self.functionXRigidity,self.functionXP_mRigidity,self.functionX_mRigidity,self.functionYPRigidity,self.functionYRigidity,self.functionYP_mRigidity,self.functionY_mRigidity,length)
		if bunch.charge() <0:
			TPB.dipoleGeneralKickStripSeperateField(bunch,firstChicaneFail,self.functionCDF,self.functionInverse,self.functionXPRigidity,self.functionXRigidity,self.functionXP_mRigidity,self.functionX_mRigidity,self.functionYPRigidity,self.functionYRigidity,self.functionYP_mRigidity,self.functionY_mRigidity,length,self.stripLength)
		elif bunch.charge()==0:
			TPB.dipoleGeneralKickStripSeperateField(bunch,secondChicaneFail,self.functionCDF,self.functionInverse,self.functionXPRigidity,self.functionXRigidity,self.functionXP_mRigidity,self.functionX_mRigidity,self.functionYPRigidity,self.functionYRigidity,self.functionYP_mRigidity,self.functionY_mRigidity,length,self.stripLength)

	