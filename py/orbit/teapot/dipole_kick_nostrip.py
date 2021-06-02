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

debugPrint=True
def printDebug(m1,m2="",m3="",m4="",m5="",m6="",m7="",m8="") :
	if debugPrint==True:
		print m1,m2,m3,m4,m5,m6,m7,m8
		
class GeneralDipoleNoStrip(BaseTEAPOT):
	def __init__(self, name = "bend no name"):
		"""
		Constructor. Creates the Dipole Combined Functions TEAPOT element .
		"""
		BaseTEAPOT.__init__(self,name)
		self.setType("child dipole teapot")
		#magnetic field in Tesla
		self.addParam("B",1.0)
		#set field direction. (0 is x-axis, pi/2 is y-xaxis)
		self.addParam("fieldDirection",1)
		#effective length of dipole
		self.addParam("effLength",.01)
		self.addParam("A1",2.47e-6)
		self.addParam("A2",4.49e9)
		
		self.functionInverse=None
		self.functionXPRigidity=None
		self.functionXRigidity=None
		self.functionXP_mRigidity=None
		self.functionX_mRigidity=None		
		
	
	def setMagneticFieldStrength(self, strength=1.0):
		#set magnetic field strength
		self.setParam("B",strength)
		
	def getMagneticFieldStrength(self):
		#returns magnetic field strength
		return self.getParam("B")
		
	def setFieldDirection(self, direction=0):
		#set the magnetic field direction
		self.setParam("fieldDirection",direction)

	def getFieldDirection(self):
		#get the magnetic field direction
		return self.getParam("fieldDirection")
		
	def setEffLength(self,effLength=0.01):
		#sets the effective length of the magnet
		self.setParam("effLength",effLength)
		
	def getEffLength(self):
		#gets the effective length of the magnet
		return self.getParam("effLength")
		
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
	def track(self, paramsDict):
		"""
		The Dipole Combined Functions TEAPOT  class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getEffLength()
		fieldDirection=self.getFieldDirection()
		bunch = paramsDict["bunch"]

		#charge <0 means first dipole
		#charge==0 means second dipole
		TPB.dipoleGeneralKickNoStrip(bunch,self.functionXPRigidity,self.functionXRigidity, length,fieldDirection)
		#if bunch.charge() <0:
			#print "hi1"
			#TPB.dipoleGeneralKickStrip(bunch,firstChicaneFail,self.functionCDF,self.functionInverse,self.functionXPRigidity,self.functionXRigidity,self.functionXP_mRigidity,self.functionX_mRigidity, length, strength,fieldDirection)
			#print "hi2"
		#elif bunch.charge()==0:
			#TPB.dipoleGeneralKickStrip(bunch,secondChicaneFail,self.functionCDF,self.functionInverse,self.functionXPRigidity,self.functionXRigidity,self.functionXP_mRigidity,self.functionX_mRigidity, length, strength,fieldDirection)

	