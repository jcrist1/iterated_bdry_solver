
mesh2 = Mesh("dom.xml")
subdomains = MeshFunction('size_t',mesh2,"dom_physical_region.xml")
boundaries = MeshFunction('size_t',mesh2,"dom_facet_region.xml")


V = FunctionSpace(mesh2, "CG",1)
Vi=vertex_to_dof_map(V)
bvals=Function(V)
pts=list()
for vertex in vertices(mesh2):
  	1+1
for vertex in vertices(mesh):
  	pts.append(vertex.point())
bbox=mesh.bounding_box_tree()
bbox.build(pts,3)
iQ=vertex_to_dof_map(Q)
for vertex in vertices(mesh2):
 	nearest_pt=bbox.compute_closest_point(vertex.point())
	value=vals.vector()[iQ[nearest_pt[0]]]
	bvals.vector()[Vi[vertex.index()]]=value[0]
File("bvals.pvd") << bvals
