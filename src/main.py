from dolfin import *
from fenics import *
from ply2fn import *
from code import *
from iterator import *
from normal_deriv import *

#comm=dolfin.mpi_comm_self()
#comm.MPI_Init_thread()

#filename="../neuron_mesh_data/neuron"
#ibdry=9
#obdry=10

filename="../head_mesh_data/head"
ibdry=10
obdry=9

#this tries to parse a ply file into a vertex mesh (is that a thing?) and spits out a dolfin xml file 

#comm=dolfin.mpi_comm_self()
#if MPI.rank(comm) == 0:
#  print(str(MPI.rank(comm))+'\n')
#gets boundary mesh and boundary values for said mesh
infile=file(filename+"_bvals.ply","r")
[mesh,Q,vals]=ply2bvals(infile)
infile.close()

mesh2 = Mesh(filename+".xml")
boundaries = MeshFunction('size_t',mesh2,filename+"_facet_region.xml")
V = FunctionSpace(mesh2, "Lagrange",1)

bvals=extend_boundary_fun(V,Q,vals)

#Compute external measurements
w=inorout(V,boundaries,bvals,ibdry,obdry,"insideout")
w.vector()[:]/=255.

nvs=get_facet_normal(V,boundaries,ibdry)
e0=Constant((1.,0.,0.))
e1=Constant((0.,1.,0.))
e2=Constant((0.,0.,1.))
nvecs=e0*nvs[0]+e1*nvs[1]+e2*nvs[2]
[delta,sol]=bval_problem_iterator(V,boundaries,w,ibdry,obdry,nvecs,65)

v=Function(V)
v.vector()[:]=w.vector()[:]-sol.vector()[:]
vsq=Function(V)
vsq.vector()[:]=v.vector()[:]*v.vector()[:]
wsq=Function(V)
wsq.vector()[:]=w.vector()[:]*w.vector()[:]
epsilon=np.sqrt(assemble(vsq*dx))
gamma=np.sqrt(assemble(wsq*dx))
error=epsilon/gamma
#File("../output/sol.pvd") << sol
File("../output/bvals.pvd") << w
#File("../output/delta.pvd") << delta
