# Nicolò Grilli
# University of Bristol
# 27 Marzo 2023

# a 2D beam element
# with force along the bar
# and can carry moments.
# It is constituted by two nodes: node1 and node2
# Notation B21 in Abaqus

import numpy as np
from element import Element
from node2d import Node2D
from material import Material

class Beam21(Element):
	
	def __init__(self,nodes,material,cross_section,moment_of_area):
		super().__init__(2,3)
		self.node1 = nodes[0]
		self.node2 = nodes[1]
		self.calculate_length()
		self.young_modulus = material.young_modulus # Young's modulus
		self.cross_section = cross_section # cross sectional area
		self.moment_of_area = moment_of_area # second moment of area
		self.calculate_stiffness_matrix()
		
	def calculate_length(self):
		self.coords1 = self.node1.coords
		self.coords2 = self.node2.coords 
		self.length = np.linalg.norm(self.coords1-self.coords2)
		if (self.length == 0):
			print("Error: zero length Truss2D2 element.")
			exit
		
	# Ke = 6x6 element stiffness matrix with rows/cols ordered as follows:
	# 1 - node 1, x DOF
	# 2 - node 1, y DOF
	# 3 - node 1, theta DOF
	# 4 - node 2, x DOF
	# 5 - node 2, y DOF
	# 6 - node 2, theta DOF
	# expressed in the global reference frame, no need for rotations
	def calculate_stiffness_matrix(self): # 6x6 stiffness matrix
		c = (self.node2.x - self.node1.x) / self.length
		s = (self.node2.y - self.node1.y) / self.length
		L = self.length
		self.stiffness_matrix = ((self.cross_section * self.young_modulus) / L) * \
								np.array([[ c*c, c*s,   0, -c*c, -c*s,  0 ],
                                          [ c*s, s*s,   0, -c*s, -s*s,  0 ],
                                          [   0,   0,   0,    0,   0,   0 ],
                                          [-c*c, -c*s,  0,  c*c,  c*s,  0 ],
                                          [-c*s, -s*s,  0,  c*s,  s*s,  0 ],
                                          [   0,   0,   0,    0,   0,   0 ]])
		self.stiffness_matrix += ((12 * self.young_modulus * self.moment_of_area) / (L*L*L)) * \
								np.array([[    s*s,  -c*s, -s*L/2,   -s*s,    c*s,  -s*L/2 ],
                                          [   -c*s,   c*c,  c*L/2,    c*s,   -c*c,   c*L/2 ],
                                          [ -s*L/2, c*L/2,  L*L/3, -s*L/2, -c*L/2,   L*L/6 ],
                                          [   -s*s,   c*s, -s*L/2,    s*s,   -c*s,   s*L/2 ],
                                          [    c*s,  -c*c, -c*L/2,   -c*s,    c*c,  -c*L/2 ],
                                          [ -s*L/2, c*L/2,  L*L/6,  s*L/2, -c*L/2,   L*L/3 ]])
