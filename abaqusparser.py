# Nicolò Grilli
# University of Bristol
# 26 Maggio 2023

# A parser for boundary conditions in Abaqus input file format

class AbaqusParser:
	
	def __init__(self,filename,BC):
		self.filename = filename
		self.BC = BC # boundary conditions object
		
	# parse boundary conditions
	def ReadBC(self):
		fid.open(self.filename,'r')
		reading_BC_flag = False
		BC_type = 'None'
		for line in fid:
			if (reading_BC_flag): # reading mode
				data = line.split()
				node_number = data[0]
				node_number = node_number.rstrip(',')
				node_number = int(node_number-1)
				dof = data[1]
				dof = dof.rstrip(',')
				dof = int(dof-1)
				magnitude = data[2]
				magnitude = magnitude.rstrip(',')
				magnitude = float(magnitude)
				
				self.BC.addDirichletBC(np.array([node_number,dof,magnitude]))

			if(line[0:9].casefold() == '*Boundary'.casefold()): # case insensitive
				reading_BC_flag = True
		
		fid.close()
