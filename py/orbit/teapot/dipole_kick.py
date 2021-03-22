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

debugPrint=True
def printDebug(m1,m2="",m3="",m4="",m5="",m6="",m7="",m8="") :
	if debugPrint==True:
		print m1,m2,m3,m4,m5,m6,m7,m8
		
class GeneralDipole(BaseTEAPOT):
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
		
	def track(self, paramsDict):
		"""
		The Dipole Combined Functions TEAPOT  class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getEffLength()
		strength=self.getMagneticFieldStrength()
		fieldDirection=self.getFieldDirection()
		bunch = paramsDict["bunch"]
		TPB.dipoleGeneralKick(bunch, length, strength,fieldDirection)
class XDipole(BaseTEAPOT):
	def __init__(self, name = "bend no name"):
		"""
		Constructor. Creates the Dipole Combined Functions TEAPOT element .
		"""
		BaseTEAPOT.__init__(self,name)
		self.setType("child dipole teapot")
		#magnetic field in Tesla
		self.addParam("B",1.0)
		#effective length of dipole
		self.addParam("effLength",.01)
		
	
	def setMagneticFieldStrength(self, strength=1.0):
		#set magnetic field strength
		self.setParam("B",strength)
		
	def getMagneticFieldStrength(self):
		#returns magnetic field strength
		return self.getParam("B")
		
	def setEffLength(self,effLength=0.01):
		#sets the effective length of the magnet
		self.setParam("effLength",effLength)
		
	def getEffLength(self):
		#gets the effective length of the magnet
		return self.getParam("effLength")
		
	def track(self, paramsDict):
		"""
		The Dipole Combined Functions TEAPOT  class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getEffLength()
		strength=self.getMagneticFieldStrength()
		bunch = paramsDict["bunch"]
		TPB.dipoleXKick(bunch, length, strength)
		
class YDipole(BaseTEAPOT):
	def __init__(self, name = "bend no name"):
		"""
		Constructor. Creates the Dipole Combined Functions TEAPOT element .
		"""
		BaseTEAPOT.__init__(self,name)
		self.setType("child dipole teapot")
		#magnetic field in Tesla
		self.addParam("B",1.0)
		#effective length of dipole
		self.addParam("effLength",.01)
		
	
	def setMagneticFieldStrength(self, strength=1.0):
		#set magnetic field strength
		self.setParam("B",strength)
		
	def getMagneticFieldStrength(self):
		#returns magnetic field strength
		return self.getParam("B")
		
	def setEffLength(self,effLength=0.01):
		#sets the effective length of the magnet
		self.setParam("effLength",effLength)
		
	def getEffLength(self):
		#gets the effective length of the magnet
		return self.getParam("effLength")
		
	def track(self, paramsDict):
		"""
		The Dipole Combined Functions TEAPOT  class implementation of
		the AccNodeBunchTracker class track(probe) method.
		"""
		length = self.getEffLength()
		strength=self.getMagneticFieldStrength()
		bunch = paramsDict["bunch"]
		TPB.dipoleYKick(bunch, length, strength)		
	