# Nicolò Grilli
# University of Bristol
# 10 Giugno 2023

# a plane strain triangle element
# Notation CPE3 used in Abaqus

import numpy as np
from constantstraintriangle import ConstantStrainTriangle

class PlaneStrainTriangle(ConstantStrainTriangle):
	
	def __init__(self,index,nodes,material):
		super().__init__(index,nodes,material)
		self.calculate_D()
		self.calculate_stiffness_matrix()

	# calculate material stiffness matrix
	def calculate_D(self):
		self.D = np.zeros(shape=(3,3))
		self.D[0,0] = 1.0 - self.poisson_ratio
		self.D[0,1] = self.poisson_ratio
		self.D[1,0] = self.poisson_ratio
		self.D[1,1] = 1.0 - self.poisson_ratio
		self.D[2,2] = 0.5 - self.poisson_ratio
		self.D *= self.young_modulus / (1.0 + self.poisson_ratio) / (1.0 - 2.0 * self.poisson_ratio)

	# calculate full element stiffness matrix
	def calculate_stiffness_matrix(self):
		DB = np.matmul(self.D, self.B)
		self.stiffness_matrix = self.area * np.matmul(np.transpose(self.B), DB)
