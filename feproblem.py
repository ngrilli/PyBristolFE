# Nicolò Grilli
# University of Bristol
# 12 Maggio 2023

# A finite element problem

import numpy as np
from node2d import Node2D
from truss2d2 import Truss2D2
from beam21 import Beam21

class FEProblem:
	
	def __init__(self,mesh,material,BC,abaqus_input_file,beam_cross_section=1,moment_of_area=1,dimensions=2):
		self.mesh = mesh
		self.material = material
		self.BC = BC # boundary conditions object
		self.abaqus_input_file = abaqus_input_file
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
		self.problem_type = self.abaqus_input_file.ReadElementType()

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
			elif (self.problem_type == 'Beam21'): # beam problem
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
			# needs to be extended for an arbitrary number of nodes per elements
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
	def apply_BC(self):
		# define force vector
		self.force_vector = np.zeros(shape=(self.dofs_per_node * self.number_of_nodes))
		for i in range(self.BC.NeumannBC.shape[0]): # cycle over Neumann BC
			self.force_vector[int(self.BC.NeumannBC[i,0] * self.dofs_per_node + self.BC.NeumannBC[i,1])] = self.BC.NeumannBC[i,2] 
		# Dirichlet BC: modify global stiffness matrix
		for i in range(self.BC.DirichletBC.shape[0]): # cycle over Dirichlet BC
			matrix_line_to_modify = int(self.BC.DirichletBC[i,0] * self.dofs_per_node + self.BC.DirichletBC[i,1])
			self.K[matrix_line_to_modify,:] = 0.0
			self.K[matrix_line_to_modify,matrix_line_to_modify] = 1.0
			self.force_vector[matrix_line_to_modify] = self.BC.DirichletBC[i,2]

	def solve(self):
		self.u = np.matmul(np.linalg.inv(self.K), self.force_vector)
