"""
Module. Includes functions that will modify the accelerator lattice by inserting the one teapot node accelerator node.
"""
# import the auxiliary classes
from orbit.utils import orbitFinalize

# import general accelerator elements and lattice
from orbit.lattice import AccLattice, AccNode, AccActionsContainer, AccNodeBunchTracker

# import Teapot dipole stripper node
from orbit.teapot import dipole_kick_nostrip
from orbit.teapot import dipole_kick_strip

# import teapot drift class and kick class
from orbit.teapot import DriftTEAPOT
from orbit.teapot import KickTEAPOT

def addDipoleStripperNode(lattice, position, stripper_node):
	"""
	It will put one Teapot collimation node in the lattice 
	"""
	length_tolerance = 0.0001
	lattice.initialize()
	position_start = position
	position_stop = position + stripper_node.getEffLength()
	(node_start_ind,node_stop_ind,z,ind) = (-1,-1, 0., 0)
	for node in lattice.getNodes():
		if(position_start >= z and position_start <= z + node.getLength()):
			node_start_ind = ind
		if(position_stop >= z and position_stop <= z + node.getLength()):
			node_stop_ind = ind
		ind += 1
		z += node.getLength()
	#-------now we check that between start and end we have only non-modified drift elements or kicks
	#-------if the space charge was added first - that is a problem. The collimation should be added first.
	for node in lattice.getNodes()[node_start_ind:node_stop_ind+1]:
		#print "debug node=",node.getName()," type=",node.getType()," L=",node.getLength()
		if(not isinstance(node,DriftTEAPOT) and not isinstance(node,KickTEAPOT)):
			print "Non-drift node=",node.getName()," type=",node.getType()," L=",node.getLength()
			orbitFinalize("We have non-drift element at the place of the foil! Stop!")
			#if(node.getNumberOfChildren() != 4):
			#print "Node=",node.getName()," type=",node.getType()," L=",node.getLength()," N children nodes=",node.getNumberOfChildren()
			#orbitFinalize("Drift element was modified with additional functionality (SC or something else)! Add collimation first! Stop!")
	# make array of nodes from foil in the center and possible two drifts if their length is more than length_tollerance [m]
	nodes_new_arr = [stripper_node,]
	drift_node_start = lattice.getNodes()[node_start_ind]
	drift_node_stop = lattice.getNodes()[node_stop_ind]	
	#check that the stripper node does not span more than 2 nodes (functionality for more not added)
	if node_stop_ind-node_start_ind>0:
		print "node_stop_ind=",node_stop_ind," node_start_ind=",node_start_ind," node_stop_ind-node_start_ind=",node_stop_ind-node_start_ind
		orbitFinalize("The stripper node spans more than 1 nodes! Stop!")
	#------now we will create two drift nodes: before the foil and after
	#------if the length of one of these additional drifts less than length_tollerance [m] we skip this drift 
	if node_stop_ind-node_start_ind==0:
		if isinstance(drift_node_start,DriftTEAPOT):
			if(position_start > lattice.getNodePositionsDict()[drift_node_start][0] +  length_tolerance):
				drift_node_start_new = DriftTEAPOT(drift_node_start.getName())
				drift_node_start_new.setLength(position_start - lattice.getNodePositionsDict()[drift_node_start][0])
				nodes_new_arr.insert(0,drift_node_start_new)
			if(position_stop < lattice.getNodePositionsDict()[drift_node_stop][1] - length_tolerance):
				drift_node_stop_new = DriftTEAPOT(drift_node_stop.getName())
				drift_node_stop_new.setLength(lattice.getNodePositionsDict()[drift_node_stop][1] - position_stop)
				nodes_new_arr.append(drift_node_stop_new)
		elif isinstance(drift_node_start,KickTEAPOT):	
			#change the strength of kicks
			totalLengthOfKicker=lattice.getNodePositionsDict()[drift_node_stop][1]-lattice.getNodePositionsDict()[drift_node_start][0]
			if(position_start > lattice.getNodePositionsDict()[drift_node_start][0] +  length_tolerance):
				drift_node_start_new = KickTEAPOT(drift_node_start.getName())
				drift_node_start_new.setLength(position_start - lattice.getNodePositionsDict()[drift_node_start][0])
				drift_node_start_new.setParam("kx",drift_node_start.getParam("kx")*drift_node_start_new.getLength()/totalLengthOfKicker)
				drift_node_start_new.setParam("ky",drift_node_start.getParam("ky")*drift_node_start_new.getLength()/totalLengthOfKicker)
				nodes_new_arr.insert(0,drift_node_start_new)
			if(position_stop < lattice.getNodePositionsDict()[drift_node_stop][1] - length_tolerance):
				drift_node_stop_new = KickTEAPOT(drift_node_stop.getName())
				drift_node_stop_new.setLength(lattice.getNodePositionsDict()[drift_node_stop][1] - position_stop)
				drift_node_stop_new.setParam("kx",drift_node_stop.getParam("kx")*drift_node_stop_new.getLength()/totalLengthOfKicker)
				drift_node_stop_new.setParam("ky",drift_node_stop.getParam("ky")*drift_node_stop_new.getLength()/totalLengthOfKicker)				
				nodes_new_arr.append(drift_node_stop_new)	
	elif node_stop_ind-node_start_ind==1:
		#this currently does work properly because of how chicane field is handled inside stripper dipole
		print "node_stop_ind-node_start_ind==1"
		orbitFinalize("node_stop_ind-node_start_ind==1! Stop!")		
		if isinstance(drift_node_start,DriftTEAPOT) and isinstance(drift_node_stop,KickTEAPOT):
			if(position_start > lattice.getNodePositionsDict()[drift_node_start][0] +  length_tolerance):
				drift_node_start_new = DriftTEAPOT(drift_node_start.getName())
				drift_node_start_new.setLength(position_start - lattice.getNodePositionsDict()[drift_node_start][0])
				nodes_new_arr.insert(0,drift_node_start_new)	
			if(position_stop < lattice.getNodePositionsDict()[drift_node_stop][1] - length_tolerance):
				drift_node_stop_new = KickTEAPOT(drift_node_stop.getName())
				drift_node_stop_new.setLength(lattice.getNodePositionsDict()[drift_node_stop][1] - position_stop)
				drift_node_stop_new.setParam("kx",drift_node_stop.getParam("kx")*drift_node_stop_new.getLength()/totalLengthOfKicker)
				drift_node_stop_new.setParam("ky",drift_node_stop.getParam("ky")*drift_node_stop_new.getLength()/totalLengthOfKicker)				
				nodes_new_arr.append(drift_node_stop_new)					
		elif isinstance(drift_node_start,KickTEAPOT) and isinstance(drift_node_stop,DriftTEAPOT):
			if(position_start > lattice.getNodePositionsDict()[drift_node_start][0] +  length_tolerance):
				drift_node_start_new = KickTEAPOT(drift_node_start.getName())
				drift_node_start_new.setLength(position_start - lattice.getNodePositionsDict()[drift_node_start][0])
				drift_node_start_new.setParam("kx",drift_node_start.getParam("kx")*drift_node_start_new.getLength()/totalLengthOfKicker)
				drift_node_start_new.setParam("ky",drift_node_start.getParam("ky")*drift_node_start_new.getLength()/totalLengthOfKicker)
				nodes_new_arr.insert(0,drift_node_start_new)
			if(position_stop < lattice.getNodePositionsDict()[drift_node_stop][1] - length_tolerance):
				drift_node_stop_new = DriftTEAPOT(drift_node_stop.getName())
				drift_node_stop_new.setLength(lattice.getNodePositionsDict()[drift_node_stop][1] - position_stop)
				nodes_new_arr.append(drift_node_stop_new)				
		else:
			print "this shouldnt be reached. isinstance(drift_node_start,DriftTEAPOT)= ",isinstance(drift_node_start,DriftTEAPOT)," isinstance(drift_node_start,KickTEAPOT)=",isinstance(drift_node_start,KickTEAPOT), " isinstance(drift_node_start,KickTEAPOT)= ",isinstance(drift_node_start,KickTEAPOT)," isinstance(drift_node_start,DriftTEAPOT)",isinstance(drift_node_start,DriftTEAPOT)
			orbitFinalize("The two nodes are not kicker and drift! Stop!")
	#------ now we will modify the lattice by replacing the found part with the new nodes
	lattice.getNodes()[node_start_ind:node_stop_ind+1] = nodes_new_arr
	# initialize the lattice
	lattice.initialize()

		


