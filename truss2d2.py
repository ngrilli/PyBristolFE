# Nicolò Grilli
# University of Bristol
# 20 Febbraio 2023

# a 2D truss element
# with force only along the bar
# This element cannot carry torque
# It is constituted by two nodes: node1 and node2
# Notation T2D2 used in Abaqus

import numpy as np
from node2d import Node2D

class Truss2D2:
	
	def __init__(self,node1,node2,cross_section,young_modulus):
		self.node1 = node1
		self.node2 = node2
		self.calculate_length()
		self.cross_section = cross_section # cross sectional area
		self.young_modulus = young_modulus # Young's modulus
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
		self.stiffness_matrix *= (self.cross_section * self.young_modulus) / self.length



