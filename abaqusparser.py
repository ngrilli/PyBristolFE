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
		for line in fid:
			if(line[0:9].casefold() == '*Boundary'.casefold()): # case insensitive
				reading_BC_flag = True
		
		fid.close()
