v=grad(u)/sqrt(inner(grad(u),grad(u))+0.1)
h=outer(v,v)
I=Identity(3)
g=I+h

psi=TrialFunction(V)
phi=TestFunction(V)
f=Constant(1.)
u0=Constant(0.)
a=dot(grad(psi),dot(grad(phi)))*dx
L=phi*f*dx
bc=DirichletBC(V,u0,boundaries,11)

psi=Function(V)
solve( a == L, psi, bc)

