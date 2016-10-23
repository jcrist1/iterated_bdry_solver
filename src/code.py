from fenics import *
# Test for PETSc
# Create mesh and define function space

def inorout(V,boundaries,bvals,ibdry,obdry,key):
  if not has_linear_algebra_backend("PETSc"):
    info("DOLFIN has not been configured with TPETSc. Exiting.")
    exit()

  parameters["linear_algebra_backend"] = "PETSc"

  #key tells us whether we solve inhomogeneous neumann on outer boundary and homogeneous dirichlet on inner for key=0, or hom neum on outer nad inhom dirichlet on inner for key=1
  # Define variational problem
  u = TrialFunction(V)
  v = TestFunction(V)
  f = Constant('0.')
  ds=Measure("ds",subdomain_data=boundaries)
  a = inner(grad(u), grad(v))*dx
  if key=='outsidein':
    #You need to thank someone for this.  The Singular Poisson Demo from Fenics
    bc=DirichletBC(V,f,boundaries,ibdry)
    L = f*v*dx - bvals*v*ds(obdry)
    # Assemble system
    A = assemble(a)
    b = assemble(L)
    bc.apply(A,b)
  
    # Solution Function
    u = Function(V)
  
    # Create Krylov solver
    solver = KrylovSolver(A, "gmres")
  
    # Create vector that spans the null space
    null_vec = Vector(u.vector())
    V.dofmap().set(null_vec, 1.0)
    null_vec *= 1.0/null_vec.norm("l2")
  
    # Create null space basis object and attach to Krylov solver
    null_space = VectorSpaceBasis([null_vec])
    solver.set_nullspace(null_space)
  
     # Orthogonalize b with respect to the null space (this gurantees that
     # a solution exists)
    null_space.orthogonalize(b);
  
  # Solve
    solver.solve(u.vector(), b)
  elif key=='insideout':
    bc=DirichletBC(V,bvals,boundaries,ibdry)
    L= f*v*dx
    u=Function(V)
    solve(a==L,u,bc)
  else:
    print("something's wrong")
  return u
