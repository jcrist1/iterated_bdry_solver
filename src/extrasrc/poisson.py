from fenics import *

mesh = Mesh("dom.xml")
subdomains = MeshFunction('size_t',mesh,"dom_physical_region.xml")
boundaries = MeshFunction('size_t',mesh,"dom_facet_region.xml")

#this tries to parse a ply file into a vertex mesh (is that a thing?) and spits out a dolfin xml file 
infile=file("bound.ply","r")
mesh2=Mesh()
editor=MeshEditor()

a=infile.readline()
while not(a.startswith("element vertex")): a=infile.readline()
tempvals=[]
b=a.split()
n=int(b[2])
while not(a.startswith("element face")): a=infile.readline()
b=a.split()
m=int(b[2])
editor.open(mesh2,2,3)
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
infile.close()
Q=FunctionSpace(mesh2,"Lagrange",1)
vals=Function(Q)
iQ=dof_to_vertex_map(Q)
for i in range(n): 
  vals.vector()[iQ[i]]=vec[i]
File("vals.pvd") << vals
# This solves poisson's equation on the inner boundary
#u=TrialFunction(Q)
#v=TestFunction(Q)
#a=inner(grad(u),grad(v))*dx
#L=vals*v*dx
#u=Function(Q)
#solve( a == L, u)
#File("vals.pvd")<< u



V = FunctionSpace(mesh, "CG",1)
Vi=vertex_to_dof_map(V)
bvals=Function(V)
pts=list()
for vertex in vertices(mesh2):
  	pts.append(vertex.point())
bbox2=mesh2.bounding_box_tree()
bbox2.build(pts,3)
for vertex in vertices(mesh):
 	nearest_pt=bbox2.compute_closest_point(vertex.point())
	value=vals.vector()[iQ[nearest_pt[0]]]
	bvals.vector()[Vi[vertex.index()]]=value[0]
File("bvals.pvd") << bvals
u = TrialFunction(V)
v = TestFunction(V)
u0= Function(V)
#vals.set_allow_extrapolation(True)
f = Constant(0.0)
a = dot(grad(u), grad(v))*dx
L = f*v*dx
bc1=DirichletBC(V,bvals,boundaries,11)
#bc2=DirichletBC(V,bvals,boundaries,13)
#bc=[bc1,bc2]
u = Function(V)
solve(a == L, u, bc1)


#File("soln.pvd") << u
