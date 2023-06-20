# Nicolò Grilli
# University of Bristol
# 14 Maggio 2023

# a generic element base class

class Element:
	
	def __init__(self,index,number_of_nodes,dofs_per_node):
		self.index = index
		self.number_of_nodes = number_of_nodes
		self.dofs_per_node = dofs_per_node
