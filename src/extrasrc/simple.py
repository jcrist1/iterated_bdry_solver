rmesh=refine(mesh,cell_markers)
rboundaries=adapt(boundaries,rmesh)
V = FunctionSpace(rmesh, "CG",1)
Vi=vertex_to_dof_map(V)
bvals=Function(V)
pts=list()
for vertex in vertices(mesh2):
  	pts.append(vertex.point())
bbox2=mesh2.bounding_box_tree()
print(vals.vector()[300])
bbox2.build(pts,3)
for vertex in vertices(rmesh):
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
bc1=DirichletBC(V,bvals,rboundaries,5)
#bc2=DirichletBC(V,bvals,boundaries,13)
#bc=[bc1,bc2]
u = Function(V)
solve(a == L, u, bc1)


