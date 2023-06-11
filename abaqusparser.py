# Nicolò Grilli
# University of Bristol
# 26 Maggio 2023

# A parser for boundary conditions in Abaqus input file format

import numpy as np

class AbaqusParser:
	
	def __init__(self,filename,material,BC):
		self.filename = filename
		self.material = material
		self.BC = BC # boundary conditions object
		
	# parse boundary conditions
	def ReadBC(self):
		fid = open(self.filename,'r')
		reading_BC_flag = False
		BC_type = 'None'
		for line in fid:
			if (line[0:1] == '*'): # find *Boundary mode
				reading_BC_flag = False
			if (reading_BC_flag): # reading mode
				data = line.split()
				node_number = data[0]
				node_number = node_number.rstrip(',')
				node_number = int(node_number)-1
				dof = data[1]
				dof = dof.rstrip(',')
				dof = int(dof)-1
				magnitude = data[2]
				magnitude = magnitude.rstrip(',')
				magnitude = float(magnitude)
				if (BC_type == 'Displacement'):
					self.BC.addDirichletBC(np.array([node_number,dof,magnitude]))
				elif (BC_type == 'Force'):
					self.BC.addNeumannBC(np.array([node_number,dof,magnitude]))
			if(line[0:9].casefold() == '*Boundary'.casefold()): # found *Boundary mode
				reading_BC_flag = True
				BC_type_string = line.split()
				BC_type_string = BC_type_string[1]
				if 'displacement' in BC_type_string.casefold():
					BC_type = 'Displacement'
				elif 'force' in BC_type_string.casefold():
					BC_type = 'Force'
				else:
					print('Error: unknown BC type parsed in input file')
					exit()
		fid.close()
		
	# parse material properties
	def ReadMaterial(self):
		fid = open(self.filename,'r')
		reading_Material_flag = False
		reading_Material_type = False
		for line in fid:
			if (reading_Material_type):
				data = line.split()
				young_modulus = data[0]
				young_modulus = young_modulus.rstrip(',')
				young_modulus = float(young_modulus)
				poisson_ratio = data[1]
				poisson_ratio = poisson_ratio.rstrip(',')
				poisson_ratio = float(poisson_ratio)
				reading_Material_flag = False
				reading_Material_type = False
			if (reading_Material_flag):
				if(line[0:8].casefold() == '*Elastic'.casefold()):
					reading_Material_type = True
				else:
					print("Error: unknown material type parsed in input file")
					exit()
			if(line[0:9].casefold() == '*Material'.casefold()): # found *Material mode
				reading_Material_flag = True
		fid.close()
		self.material.young_modulus = young_modulus
		self.material.poisson_ratio = poisson_ratio

	# parse element type
	def ReadElementType(self):
		fid = open(self.filename,'r')
		for line in fid:
			if (line[0:8].casefold() == '*Element'.casefold()):
				data = line.split()
				data = data[1]
				data = data.lstrip('type=')
				if (data == 'T2D2'):
					problem_type = 'Truss2D2'
				elif (data == 'B21'):
					problem_type = 'Beam21'
				elif (data == 'S3'):
					problem_type = 'PlainStrainTriangle'
				elif (data == 'CPS3'):
					problem_type = 'PlainStressTriangle'
				else:
					print("Error: unknown element type parsed in input file")
					exit()
				break
		fid.close()
		return problem_type
