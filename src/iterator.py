from dolfin import *
import numpy as np
from code import *
#N is the number of iterations
def bval_problem_iterator(V,boundaries,w,ibdry,obdry,nvec,N):
  solfile=File("../output/sol.pvd")
  deltafile=File("../output/delta.pvd")
  ds=Measure("ds",subdomain_data=boundaries)
  TOL=0.00000000001
  wsq=float(assemble(w*w*dx))
  wsqbd=float(assemble(w*w*ds(ibdry)))
  outfile=open('../output/error.txt','w')
  outfile.write('iter,\t relerr\n')
  log2=np.log(2.)
  bvals=Function(V)
  bvals.vector()[:]=w.vector()[:]
  f=Function(V)
  f.vector()[:]=1.
  ds=Measure("ds",subdomain_data=boundaries)
  area=assemble(f*ds(obdry))
  vol=assemble(f*dx)
  sol=Function(V)
  avg=float(assemble(bvals*ds(obdry))/area)
  bvals.vector()[:]-=avg
  sol.vector()[:]+=avg
  for i in range(N):
    w1=inorout(V,boundaries,bvals,ibdry,obdry,'outsidein')
    u1=project(dot(grad(w1),nvec),V)
    w2=inorout(V,boundaries,u1,ibdry,obdry,'insideout')
    #w3=inorout(V,boundaries,w2,ibdry,obdry,'outsidein')
    c=float(assemble(w2*bvals*ds(obdry)))
    print(c)
    #+float(assemble(inner(grad(w1),grad(w3))*dx))
    bsq=float(assemble(w2*w2*ds(obdry)))
    print(bsq)
    #+float(assemble(inner(grad(w3),grad(w3))*dx))
    tau=c/bsq
    sol.vector()[:] += tau*w2.vector()[:]
    bvals.vector()[:]-=tau*w2.vector()[:]
    avg=float(assemble(bvals*ds(obdry))/area)
    sol.vector()[:]+=avg
    bvals.vector()[:]-=avg
    log_2i=np.log(float(i))/log2
    if abs(log_2i-round(log_2i))<TOL:
       delta=Function(V)
       delta.vector()[:]=w.vector()[:]-sol.vector()[:]
       delta.rename('error','error')
       deltasq=float(assemble(delta*delta*dx))
       deltasqbd=float(assemble(delta*delta*ds(ibdry)))
       error=np.sqrt(deltasq/wsq)
       errorbd=np.sqrt(deltasqbd/wsqbd)
       sol.rename('solution','solution')
       solfile <<sol
       deltafile<<delta
       outfile.write(str(i)+',\t'+str(error)+',\t'+str(errorbd)+'\n')
  outfile.close()
  return [bvals,sol] 
