# Nicolò Grilli
# University of Bristol
# 2 Giugno 2023

# VTK file writer class

import meshio

class VTKwriter:
	
	def _init_(self,filename,mesh,feproblem):
		self.filename = filename # output file name
		self.mesh = mesh
		self.feproblem = feproblem
		
	# create output mesh with displacement nodal field
	def load_output_field(self):
		return 1

	def write_output(self):
		self.output_mesh.write(self.filename,file_format="vtk")
		
