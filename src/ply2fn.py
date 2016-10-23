from fenics import *
from dolfin import *

#gets boundary mesh and boundary values for said mesh

def ply2bvals(infile):
  mesh=Mesh()
  editor=MeshEditor()
  a=infile.readline()
  while not(a.startswith("element vertex")): a=infile.readline()
  tempvals=[]
  b=a.split()
  n=int(b[2])
  while not(a.startswith("element face")): a=infile.readline()
  b=a.split()
  m=int(b[2])
  editor.open(mesh,2,3)
  editor.init_vertices(n)
  editor.init_cells(m)
  
  while a!="end_header\n" and a!="" : a=infile.readline()
  if a=="": print("help EOF\n")
  a=infile.readline()
  b=a.split()
  vec=[]
  for i in range(n):
  	editor.add_vertex(i,float(b[0]),float(b[1]),float(b[2]))
  	vec=vec+[float(b[8])]
  	a=infile.readline()
  	b=a.split(None)
  for i in range(m):
    	editor.add_cell(i,int(b[1]),int(b[2]),int(b[3]))
  	a=infile.readline()
  	b=a.split()
  editor.close()
  Q=FunctionSpace(mesh,"Lagrange",1)
  vals=Function(Q)
  iQ=vertex_to_dof_map(Q)
  for i in range(n): 
    vals.vector()[iQ[i]]=vec[i]
  return [mesh,Q,vals]

#V = FunctionSpace(mesh2, "CG",1)
def extend_boundary_fun(V,Q,vals):
  mesh=Q.mesh()
  mesh2=V.mesh()
  Vi=vertex_to_dof_map(V)
  bvals=Function(V)
  pts=list()
  for vertex in vertices(mesh):
    	pts.append(vertex.point())
  bbox=mesh.bounding_box_tree()
  bbox.build(pts,3)
  iQ=vertex_to_dof_map(Q)
  for vertex in vertices(mesh2):
   	nearest_pt=bbox.compute_closest_point(vertex.point())
  	value=vals.vector()[iQ[nearest_pt[0]]]
  	bvals.vector()[Vi[vertex.index()]]=value[0]
  return bvals
