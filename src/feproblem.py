# Nicolò Grilli
# University of Bristol
# 12 Maggio 2023

# A finite element problem

import numpy as np
from node2d import Node2D
from truss2d2 import Truss2D2
from beam21 import Beam21
from planestresstriangle import PlaneStressTriangle
from planestraintriangle import PlaneStrainTriangle

class FEProblem:
	
	def __init__(self,mesh,material,BC,abaqus_input_file):
		self.mesh = mesh
		self.material = material
		self.BC = BC # boundary conditions object
		self.abaqus_input_file = abaqus_input_file
		
	def calculate_number_of_nodes(self):
		self.number_of_nodes = self.mesh.points.__len__()
		self.nodes = []
		# create a list of node objects
		for n in range(self.number_of_nodes):
			self.nodes.append(Node2D(n,self.mesh.points[n][:]))
		
	def find_element_type(self):
		self.element_type = list(self.mesh.cells_dict.keys())[0]
		self.problem_type = self.abaqus_input_file.problem_type

	def calculate_number_of_elements(self):
		self.number_of_elements = len(self.mesh.cells_dict[self.element_type])
		self.elements = []
		# create a list of element objects
		for elem in range(self.number_of_elements):
			nodes = [] # list of nodes of this element
			if (self.problem_type == 'Truss2D2'): # pin-jointed bar problem
				for i in range(2):
					nodes.append(self.nodes[self.mesh.cells_dict['line'][elem][i]])
				self.elements.append(Truss2D2(elem,nodes,self.material))
			elif (self.problem_type == 'Beam21'): # beam problem
				for i in range(2):
					nodes.append(self.nodes[self.mesh.cells_dict['line'][elem][i]])
				self.elements.append(Beam21(elem,nodes,self.material))
			elif (self.problem_type == 'PlaneStrainTriangle'): # plane strain triangle
				for i in range(3):
					nodes.append(self.nodes[self.mesh.cells_dict['triangle'][elem][i]])
				self.elements.append(PlaneStrainTriangle(elem,nodes,self.material))
			elif (self.problem_type == 'PlaneStressTriangle'): # plane stress triangle
				for i in range(3):
					nodes.append(self.nodes[self.mesh.cells_dict['triangle'][elem][i]])
				self.elements.append(PlaneStressTriangle(elem,nodes,self.material))
				
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
			# cycle over pairs of nodes in this element and their degrees of freedom
			for n1 in range(len(elem.nodes)):
				for n2 in range(len(elem.nodes)):
					for i in range(elem.dofs_per_node):
						for j in range(elem.dofs_per_node):
							# find global matrix indices ii and jj for this Ke entry
							node1 = elem.nodes[n1]
							node2 = elem.nodes[n2]
							ii = node1.index * elem.dofs_per_node + i
							jj = node2.index * elem.dofs_per_node + j
							self.K[ii,jj] += Ke[i + n1 * elem.dofs_per_node, j + n2 * elem.dofs_per_node]

	# apply boundary conditions
	def apply_BC(self):
		# define force vector
		self.force_vector = np.zeros(shape=(self.dofs_per_node * self.number_of_nodes))
		for i in range(self.BC.NeumannBC.shape[0]): # cycle over Neumann BC
			self.force_vector[int(self.BC.NeumannBC[i,0] * self.dofs_per_node + self.BC.NeumannBC[i,1])] = self.BC.NeumannBC[i,2] 
		# store original global stiffness matrix for force output
		self.K_original = np.copy(self.K)
		# Dirichlet BC: modify global stiffness matrix
		for i in range(self.BC.DirichletBC.shape[0]): # cycle over Dirichlet BC
			matrix_line_to_modify = int(self.BC.DirichletBC[i,0] * self.dofs_per_node + self.BC.DirichletBC[i,1])
			self.K[matrix_line_to_modify,:] = 0.0
			self.K[matrix_line_to_modify,matrix_line_to_modify] = 1.0
			self.force_vector[matrix_line_to_modify] = self.BC.DirichletBC[i,2]

	def solve(self):
		self.u = np.matmul(np.linalg.inv(self.K), self.force_vector)
		self.force_vector_output = np.matmul(self.K_original, self.u)

	# calculate strain for triangle elements
	def calculate_strain(self):
		self.Exx = []
		self.Eyy = []
		self.Exy = []
		self.Ezz = []
		if (self.problem_type == 'PlaneStrainTriangle' or self.problem_type == 'PlaneStressTriangle'):
			for elem in self.elements:
				u_elem_index = 0 # index of displacement solution in this element
				u_elem = np.zeros(shape=(elem.number_of_nodes * elem.dofs_per_node))
				for n in range(elem.number_of_nodes):
					node = elem.nodes[n]
					for i in range(elem.dofs_per_node):
						ii = node.index * elem.dofs_per_node + i
						# displacement solution in this elements
						u_elem[u_elem_index] = self.u[ii]
						u_elem_index += 1
				strain_elem = np.matmul(elem.B, u_elem) # strain vector in this element
				self.Exx.append(strain_elem[0])
				self.Eyy.append(strain_elem[1])
				self.Exy.append(strain_elem[2])
				if (self.problem_type == 'PlaneStrainTriangle'):
					self.Ezz.append(0.0)
				if (self.problem_type == 'PlaneStressTriangle'):
					nu = self.material.poisson_ratio
					self.Ezz.append(-nu / (1.0 - nu) * (strain_elem[0] + strain_elem[1]))
				
	# calculate stress for triangle elements
	def calculate_stress(self):
		self.Sxx = []
		self.Syy = []
		self.Sxy = []
		self.Szz = []
		if (self.problem_type == 'PlaneStrainTriangle' or self.problem_type == 'PlaneStressTriangle'):
			for elem in self.elements:
				strain_elem = np.array([self.Exx[elem.index], self.Eyy[elem.index], self.Exy[elem.index]])
				stress_elem = np.matmul(elem.D, strain_elem)
				self.Sxx.append(stress_elem[0])
				self.Syy.append(stress_elem[1])
				self.Sxy.append(stress_elem[2])
				if (self.problem_type == 'PlaneStressTriangle'):
					self.Szz.append(0.0)
				if (self.problem_type == 'PlaneStrainTriangle'):
					E = self.material.young_modulus
					nu = self.material.poisson_ratio
					self.Szz.append(E * nu / (1.0 + nu) / (1.0 - 2.0 * nu) * (strain_elem[0] + strain_elem[1]))
		
