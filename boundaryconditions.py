# Nicolò Grilli
# University of Bristol
# 26 Maggio 2023

# Boundary conditions

import numpy as np

class BoundaryConditions:
	
	# data structure is a matrix with 3 columns for each line, the 3 columns mean:
	# node number, degree of freedom, magnitude
	# zero index notation is used
	def __init__(self):
		self.DirichletBC = np.zeros(shape=(0,3))
		self.NeumannBC = np.zeros(shape=(0,3))
		
	def addDirichletBC(self,newBC):
		self.DirichletBC = np.vstack(self.DirichletBC, newBC)

	def addNeumannBC(self,newBC):
		self.NeumannBC = np.vstack(self.NeumannBC, newBC)
