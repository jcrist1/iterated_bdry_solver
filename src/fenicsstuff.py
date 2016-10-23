from dolfin import *
mesh = Mesh("dom.xml")
subdomains = MeshFunction('size_t',mesh,"dom_physical_region.xml")
boundaries = MeshFunction('size_t',mesh,"dom_facet_region.xml")
p=Point()
mesh.closest_point(p)
#V = FunctionSpace(mesh, "CG",1)
#u = TrialFunction(V)
#v = TestFunction(V)
#f = Constant(3.0)
#a = dot(grad(u), grad(v))*dx
#L = f*v*dx
#u0=Constant(0.0)
#u1=Constant(4.0)
#bc1=DirichletBC(V,u0,boundaries,8)
#bc2=DirichletBC(V,u1,boundaries,9)
#bc=[bc1,bc2]
#u = Function(V)
#solve(a == L, u, bc)
##plot(subdomains,interactive=True)
##plot(mesh,interactive=True)
##plot(boundaries,interactive=True)
#plot(u, interactive=True)
#
#
#File("soln.pvd") << u
