from dolfin import  *
import numpy as np

def get_facet_normal(V,boundaries,bval):
  mesh=V.mesh()
  mesh.init(0,2)
  mesh.init(2,3)
  normal0=Function(V)
  normal1=Function(V)
  normal2=Function(V)
  normal=[normal0,normal1,normal2]
  coords=mesh.coordinates()
  cells=mesh.cells()
  iV=vertex_to_dof_map(V)
  facetlist=[] 
  facetcells=[]
  for fcts in facets(mesh):
     facetlist+=[fcts.entities(0)]
     facetcells+=[fcts.entities(3)]
  nvs=[]
  for vrtx in vertices(mesh):
    key=0
    nv=[0.,0.,0.]
    for fctnum in vrtx.entities(2):
      #print(fctnum)
      fctpts=facetlist[fctnum]
      #print(fctpts)
      #print(boundaries[int(fctnum)])
      if boundaries[int(fctnum)]==bval:
	if key==0:
	  key=1
	cllpts=cells[facetcells[fctnum]].flatten()
        j=0
	#print(cllpts)
	#print(fctpts)
	for i in (0,1,2,3):
	  if cllpts[i] not in fctpts:
	    j=i
        vec1=coords[fctpts[1]]-coords[fctpts[0]]
	vec2=coords[fctpts[2]]-coords[fctpts[0]]
	vec3=coords[cllpts[j]]-coords[fctpts[0]]
	n=np.cross(vec1,vec2)
	if np.dot(n,vec3)>=0:
	  n=-n
        nv+=n
    if key==1:
      normalnorm =np.sqrt(np.dot(nv,nv))
      nv /=normalnorm
      nvs=nvs+[nv]
    for k in (0,1,2):
      normal[k].vector()[iV[vrtx.index()]]=nv[k]
  return normal
#nvs=get_facet_normal(V,boundaries,ibdry)
#e0=Constant((1.,0.,0.))
#e1=Constant((0.,1.,0.))
#e2=Constant((0.,0.,1.))
#nvec=nvs[0]*e0+nvs[1]*e1+nvs[2]*e2
#

