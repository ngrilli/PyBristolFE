# Nicolò Grilli
# University of Bristol
# 12 Maggio 2023

# A finite element problem

import numpy as np
from material import Material
from node2d import Node2D
from truss2d2 import Truss2D2
from beam21 import Beam21

class FEProblem:
	
	def __init__(self,problem_type,mesh,material,BC,beam_cross_section=1,moment_of_area=1,dimensions=2):
		self.problem_type = problem_type
		self.mesh = mesh
		self.material = material
		self.beam_cross_section = beam_cross_section
		self.moment_of_area = moment_of_area
		self.dimensions = dimensions
		
	def calculate_number_of_nodes(self):
		self.number_of_nodes = self.mesh.points.__len__()
		self.nodes = []
		# create a list of node objects
		for n in range(self.number_of_nodes):
			self.nodes.append(Node2D(n,self.mesh.points[n][:]))
		
	def find_element_type(self):
		self.element_type = list(self.mesh.cells_dict.keys())[0]

	def calculate_number_of_elements(self):
		self.number_of_elements = len(self.mesh.cells_dict[self.element_type])
		# How to get .inp element type with meshio?
		self.elements = []
		# create a list of element objects
		for elem in range(self.number_of_elements):
			node1 = self.nodes[self.mesh.cells_dict['line'][elem][0]]
			node2 = self.nodes[self.mesh.cells_dict['line'][elem][1]]
			if (self.problem_type == 'Truss2D2'): # pin-jointed bar problem
				self.elements.append(Truss2D2([node1,node2],self.material,self.beam_cross_section))
			if (self.problem_type == 'Beam21'): # beam problem
				self.elements.append(Beam21([node1,node2],self.material,self.beam_cross_section,self.moment_of_area))
				
	# would be better to sum over elements in case of mesh with different element types 
	def calculate_number_of_dofs_per_node(self):
		first_element = self.elements[0]
		self.dofs_per_node = first_element.dofs_per_node

	# assemble global stiffness matrix
	def calculate_global_stiffness_matrix(self):
		K_dimensions = self.dofs_per_node * self.number_of_nodes 
		self.K = np.zeros(shape=(K_dimensions,K_dimensions))
		for elem in self.elements:
			Ke = elem.stiffness_matrix
			Ke_dimensions = len(Ke)
			node1 = elem.node1
			node2 = elem.node2
			# find global matrix indices ii and jj for each element
			for i in range(elem.dofs_per_node):
				for j in range(elem.dofs_per_node):
					# (node1 , node1) submatrix
					ii = node1.index * elem.dofs_per_node + i
					jj = node1.index * elem.dofs_per_node + j
					self.K[ii,jj] += Ke[i,j]
					# (node1 , node2) submatrix
					ii = node1.index * elem.dofs_per_node + i
					jj = node2.index * elem.dofs_per_node + j
					self.K[ii,jj] += Ke[i,j+elem.dofs_per_node]
					# (node2 , node1) submatrix
					ii = node2.index * elem.dofs_per_node + i
					jj = node1.index * elem.dofs_per_node + j
					self.K[ii,jj] += Ke[i+elem.dofs_per_node,j]
					# (node2 , node2) submatrix
					ii = node2.index * elem.dofs_per_node + i
					jj = node2.index * elem.dofs_per_node + j
					self.K[ii,jj] += Ke[i+elem.dofs_per_node,j+elem.dofs_per_node]

	# apply boundary conditions
	# a BC is a matrix with 8 or 12 columns for Truss2D2 or Beam21
	# 1st column = value of force_x on node
	# 2nd column = is force_x defined on this node (1 = True, 0 = False)
	# 3rd column = value of force_y on node
	# 4th column = is force_y defined on this node (1 = True, 0 = False)
	# 5th column = value of displacement_x on node
	# 6th column = is displacement_x defined on this node (1 = True, 0 = False)
	# 7th column = value of displacement_y on node
	# 8th column = is displacement_y defined on this node (1 = True, 0 = False)
	# 9th column = value of torque defined on node
	# 10th column = is torque defined on this nodes (1 = True, 0 = False)
	# 11th column = value of rotation defined on node
	# 12th column = is rotation defined on this node (1 = True, 0 = False)
	def apply_BC(self):
		# split stiffness matrix into known/unknown submatrices
		return 1
		
	def solve(self):
		return 1
