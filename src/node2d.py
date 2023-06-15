# Nicolò Grilli
# University of Bristol
# 20 Febbraio 2023

# A 2D node

class Node2D:

	def __init__(self,index,coords):
		self.index = int(index)
		self.coords = coords
		if (len(coords) != 2):
			print("Warning: wrong array length in Node2D constructor.")
		self.x = coords[0]
		self.y = coords[1]

