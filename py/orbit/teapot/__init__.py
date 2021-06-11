## \namespace orbit::teapot
## \brief Python classes for TEAPOT elements.
##
## These classes use teapot_base C++ wrappers

from teapot import TEAPOT_Lattice
from teapot import BaseTEAPOT
from teapot import BendTEAPOT
from teapot import DriftTEAPOT
from teapot import FringeFieldTEAPOT
from teapot import KickTEAPOT
from teapot import MultipoleTEAPOT
from teapot import QuadTEAPOT
from teapot import RingRFTEAPOT
from teapot import SolenoidTEAPOT
from teapot import TiltTEAPOT


from dipole_kick import GeneralDipole
from dipole_kick import XDipole
from dipole_kick import YDipole

from dipole_kick_strip import GeneralDipoleStrip
from dipole_kick_nostrip import GeneralDipoleNoStrip

from dipole_kick_strip_separate_field_components import GeneralDipoleStripSeperateField
from dipole_kick_nostrip_separate_field_components import GeneralDipoleNoStripSeperateField

from teapot import TPB
from dipoleStripperLatticeModifications import addDipoleStripperNode

from teapot_matrix_lattice import TEAPOT_MATRIX_Lattice

__all__ = []
__all__.append("TEAPOT_Lattice")
__all__.append("BaseTEAPOT")
__all__.append("DriftTEAPOT")
__all__.append("BunchWrapTEAPOT")
__all__.append("BendTEAPOT")
__all__.append("QuadTEAPOT")
__all__.append("MultipoleTEAPOT")
__all__.append("SolenoidTEAPOT")
__all__.append("KickTEAPOT")
__all__.append("RingRFTEAPOT")
__all__.append("FringeFieldTEAPOT")
__all__.append("TiltTEAPOT")
__all__.append("TPB")
__all__.append("TEAPOT_MATRIX_Lattice")
__all__.append("GeneralDipole")
__all__.append("XDipole")
__all__.append("YDipole")

