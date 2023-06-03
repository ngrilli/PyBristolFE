# Nicolò Grilli
# University of Bristol
# 2 Giugno 2023

# VTK file writer class

import meshio

class VTKwriter:
	
	def __init__(self,filename,mesh,feproblem):
		self.filename = filename # output file name
		self.mesh = mesh
		self.feproblem = feproblem
		
	# create output mesh with displacement nodal field
	def load_output_field(self):
		solution_dict = {} # solution dictionary for output
		for dof in range(self.feproblem.dofs_per_node):
			solution_array = []
			for node in range(self.feproblem.number_of_nodes):
				solution_array.append(self.feproblem.u[node*self.feproblem.dofs_per_node + dof])
			if (dof == 0):
				solution_dict['ux'] = solution_array
			elif (dof == 1):
				solution_dict['uy'] = solution_array
			elif (dof == 2):
				solution_dict['utheta'] = solution_array
		self.output_mesh = meshio.Mesh(self.mesh.points,self.mesh.cells,point_data=solution_dict,)

	def write_output(self):
		self.output_mesh.write(self.filename,file_format="vtk")
		
