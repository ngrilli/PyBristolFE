# Nicolò Grilli
# University of Bristol
# 6 Giugno 2023

# constant strain triangle
# generic element base class

from element import Element
import numpy as np

class ConstantStrainTriangle(Element):
	
	def __init__(self,nodes):
		super().__init__(3,2)
		self.node1 = nodes[0]
		self.node2 = nodes[1]
		self.node3 = nodes[2]
		self.calculate_area()
		self.calculate_dshape_function()
		
	# calculate element area
	def calculate_area(self):
		self.x1 = self.node1.coords[0]
		self.x2 = self.node2.coords[0]
		self.x3 = self.node3.coords[0]
		self.y1 = self.node1.coords[1]
		self.y2 = self.node2.coords[1]
		self.y3 = self.node3.coords[1]
		self.z1 = self.node1.coords[2]
		self.z2 = self.node2.coords[2]
		self.z3 = self.node3.coords[2]
		self.area = 0.5 * (self.x2 * self.y3 - self.x3 * self.y2 + \
		                   self.x3 * self.y1 - self.x1 * self.y3 + \
		                   self.x1 * self.y2 - self.x2 * self.y1)

	# calculate B matrix
	def calculate_dshape_function(self):
		self.B = np.zeros(shape=(self.number_of_nodes,self.number_of_nodes * self.dofs_per_node))
		self.B[0,0] = self.y2 - self.y3
		self.B[0,2] = self.y3 - self.y1
		self.B[0,4] = self.y1 - self.y2
		self.B[1,1] = self.x3 - self.x2
		self.B[1,3] = self.x1 - self.x3
		self.B[1,5] = self.x2 - self.x1
		self.B[2,0] = self.x3 - self.x2
		self.B[2,1] = self.y2 - self.y3
		self.B[2,2] = self.x1 - self.x3
		self.B[2,3] = self.y3 - self.y1
		self.B[2,4] = self.x2 - self.x1
		self.B[2,5] = self.y1 - self.y2
		# divide by twice the area
		self.B = self.B / (2.0 * self.area)
