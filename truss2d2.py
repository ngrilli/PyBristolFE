# Nicolò Grilli
# University of Bristol
# 20 Febbraio 2023

# a 2D truss element
# with force only along the bar
# This element cannot carry torque
# It is constituted by two nodes: node1 and node2
# Notation T2D2 used in Abaqus
# Also called pin-jointed bar

import numpy as np
from element import Element
from node2d import Node2D
from material import Material

class Truss2D2(Element):
	
	def __init__(self,nodes,material):
		super().__init__(2,2)
		self.nodes = nodes
		self.node1 = nodes[0]
		self.node2 = nodes[1]
		self.calculate_length()
		self.young_modulus = material.young_modulus # Young's modulus
		self.beam_cross_section = material.beam_cross_section # cross sectional area
		self.calculate_stiffness_matrix()
		
	def calculate_length(self):
		self.coords1 = self.node1.coords
		self.coords2 = self.node2.coords 
		self.length = np.linalg.norm(self.coords1-self.coords2)
		if (self.length == 0):
			print("Error: zero length Truss2D2 element.")
			exit
		
	# Ke = 4x4 element stiffness matrix with rows/cols ordered as follows:
	# 1 - node 1, x DOF
	# 2 - node 1, y DOF
	# 3 - node 2, x DOF
	# 4 - node 2, y DOF
	def calculate_stiffness_matrix(self): # 4x4 stiffness matrix
		c = (self.node2.x - self.node1.x) / self.length
		s = (self.node2.y - self.node1.y) / self.length
		self.stiffness_matrix = np.array([[ c*c,  c*s, -c*c, -c*s],
                                          [ c*s,  s*s, -c*s, -s*s],
                                          [-c*c, -c*s,  c*c,  c*s],
                                          [-c*s, -s*s,  c*s,  s*s]])
		self.stiffness_matrix *= (self.beam_cross_section * self.young_modulus) / self.length



