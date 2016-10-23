from dolfin import *

#Load the Mesh

#Get the boundary values

#Compute the forward solution

#Iteration 
	#compute adjoint squared



#def forward solver has functionspace property, and one method:solve which takes interior bvals arg.  dirichlet bvals on interior and 0 neumann on exterior.

#def backwards solver has function space property, boundaries (MeshFunction) property, interior boundary, and exterior property, and one method: solve which takes exterior neumann values, applying 0 dirichlet values on interior

#still need to get neumann data out somehow.
