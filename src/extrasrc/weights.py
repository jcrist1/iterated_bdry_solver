#this tries to parse a ply file into a vertex mesh (is that a thing?) and spits out a dolfin xml file 
from dolfin import *
infile=file("idomain.ply","r")
outfile1=File("weight_vertices.xml")
outfile2=File("weights.xml")
mesh=Mesh()
editor=MeshEditor()

a=infile.readline()
while not(a.startswith("element vertex")): a=infile.readline()
tempvals=[]
b=a.split()
n=int(b[2])
editor.open(mesh,0,3)
editor.init_vertices(n)
editor.init_cells(n)
vals=MeshFunctionDouble(mesh)
vals.init(0)
while a!="end_header\n" and a!="" : a=infile.readline()
if a=="": print("help EOF\n")
a=infile.readline()
b=a.split()
for i in range(n):
	editor.add_vertex(i,float(b[0]),float(b[1]),float(b[2]))
	vals.set_value(i,float(b[10]))
	a=infile.readline()
	b=a.split(None)
editor.close()
infile.close()
outfile1 << mesh
outfile2 << vals
plot(vals, interactive=True)



#<?xml version="1.0" encoding="UTF-8"?>
#
#<dolfin xmlns:dolfin="http://www.fenicsproject.org">
#  <mesh celltype="vertex" dim="0">
#    <vertices size="8178">
#      <vertex index="0" x="8.8531053066253662e-01" y="-6.2765264511108398e-01" z="1.8571165204048160e-01"/>
#      <vertex index="1" x="8.4808081388473511e-01" y="-6.7841160297393799e-01" z="1.6500797867774961e-01"/>
#  </mesh>
#</dolfin>
#x	y		z	nx	ny	nz	s	 t	  r   g   b
#0.885311 -0.627653 0.185712 0.628498 -0.576556 0.522050 0.378894 0.966613 255 255 255

