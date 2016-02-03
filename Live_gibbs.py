from dolfin import *

N = 12


class F(Expression):
	def eval(selv,v,x):
		if x[0]<0.66666 and x[0]>0.3333: v[0] = 1
		else: v[0] = 0

Ns = [10,20,40,80,160]
for N in Ns:


	mesh = UnitIntervalMesh(N)

	f = F()
	V = FunctionSpace(mesh, "Lagrange", 1)
	V3 = FunctionSpace(mesh, "Lagrange", 1+3)
	f3 = project(f,V3)
	u = TrialFunction(V)
	v = TestFunction(V)

	a = u*v*dx
	L = f*v*dx
	U = Function(V)
	solve(a==L,U)
	#plot(U)
	#interactive()
	L2_error = sqrt(assemble((U-f3)**2*dx))
	H1_error = sqrt(assemble(inner(grad(U-f3),grad(U-f3))*dx + (U-f3)**2*dx) )
	print "N ", N, "L2_error ",L2_error
	print "N ", N ,"H1_error ", H1_error
