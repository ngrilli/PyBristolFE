# Nicolò Grilli
# University of Bristol
# 26 Maggio 2023

# A parser for boundary conditions in Abaqus input file format

import numpy as np

class AbaqusParser:
	
	def __init__(self,filename,BC):
		self.filename = filename
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
		return 1
