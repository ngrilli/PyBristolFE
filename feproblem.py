# Nicolò Grilli
# University of Bristol
# 12 Maggio 2023

# A finite element problem

import numpy as np

class FEProblem:
	
	def __init__(self,problem_type,mesh,young_modulus,beam_cross_section=1,moment_of_area=1,dimensions=2):
		self.problem_type = problem_type
		self.mesh = mesh
		self.young_modulus = young_modulus
		self.beam_cross_section = beam_cross_section
		self.moment_of_area = moment_of_area
		self.dimensions = dimensions
		
	def calculate_number_of_nodes():
		self.number_of_nodes = mesh.points.__len__()
		self.nodes = []
		# create a list of node objects
		for n in range(self.number_of_nodes):
			self.nodes.append(Node2D(n,mesh.points[n][:]))
		
	def find_element_type():
		self.element_type = list(mesh.cells_dict.keys())[0]

	def calculate_number_of_elements():
		self.number_of_elements = len(mesh.cells_dict[self.element_type])
		# How to get .inp element type with meshio?
		self.elements = []
		# create a list of element objects
		for elem in range(self.number_of_elements):
			node1 = self.nodes[mesh.cells_dict['line'][elem][0]]
			node2 = self.nodes[mesh.cells_dict['line'][elem][1]]
			if (self.problem_type == 'Truss2D2'): # pin-jointed bar problem
				self.elements.append(Truss2D2(node1,node2,self.beam_cross_section,self.young_modulus))
			if (self.problem_type == 'Beam21'): # beam problem
				self.elements.append(Beam21(node1,node2,self.beam_cross_section,self.young_modulus,self.moment_of_area))

	# assemble global stiffness matrix
	def calculate_global_stiffness_matrix():
		K_dimensions = self.dimensions * self.number_of_nodes 
		self.K = np.zeros(shape=(K_dimensions,K_dimensions))
		for elem in self.elements:
			Ke = elem.stiffness_matrix
			Ke_dimensions = len(Ke)
			# find global matrix indices for each element
			# DA CONTROLLARE -> altri due loop sui nodi necessari
			node1 = elem.node1
			node2 = elem.node2
			for i in range(elem.dofs_per_node):
				for j in range(elem.dofs_per_node):
					# (node 1 , node 1) submatrix
					ii = node1.index * elem.dofs_per_node + i
					jj = node1.index * elem.dofs_per_node + j
					self.K[ii,jj] += Ke[i,j]
					# second node
					j = node2.index * elem.dofs_per_node + dof


	# apply boundary conditions
	def apply_BC():
		return 1
