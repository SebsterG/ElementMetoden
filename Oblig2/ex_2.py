from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

set_log_active(False)
def solve_me(N ,mu,lamb):
	mesh = UnitSquareMesh(N,N)
	V = VectorFunctionSpace(mesh,"CG",3)
	W = VectorFunctionSpace(mesh,"CG",4)
	Q = FunctionSpace(mesh,"CG",2)
	VQ = V*Q
	up = TrialFunction(VQ)
	u, p = split(up)
	vq = TestFunction(VQ)
	v, q = split(vq)

	class Boundary(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary
	u_e = Expression(("pi*x[0]*cos(pi*x[0]*x[1])","-pi*x[1]*cos(pi*x[0]*x[1])"))
	u_exact = project(u_e,W)
	#plot(u_e_plot,interactive=True)

	f = Expression(("pi*pi*(pi*x[0]*(x[0]*x[0]+x[1]*x[1])*cos(pi*x[0]*x[1])+2*x[1]*sin(pi*x[0]*x[1]))",\
	                 " -pi*pi*(pi*x[1]*(x[0]*x[0]+x[1]*x[1])*cos(pi*x[0]*x[1])+2*x[0]*sin(pi*x[0]*x[1]))"))
	boundary = Boundary()
	bound = FacetFunction("size_t",mesh)
	boundary.mark(bound,1)
	#plot(bound,interactive=True)

	bc1 = DirichletBC(VQ.sub(0), u_e, boundary)
	#bc2 = DirichletBC(Q, 0, boundary)
	bcs = [bc1]
	mu = Constant(mu); lamb = Constant(lamb)

	a1 = mu*inner(grad(u),grad(v))*dx + p*div(v)*dx
	a2 = - (1.0/lamb)*p*q*dx  + dot(div(u),q)*dx
	L = inner(f,v)*dx
	u_p_ = Function(VQ)
	solve(a1+a2==L,u_p_,bcs)
	u_, p_ = u_p_.split()
	#plot(u_,interactive=True)
	return u_, u_exact,mesh

def solve_me_bad(N,mu,lamb):
	mesh = UnitSquareMesh(N,N)
	V = VectorFunctionSpace(mesh,"CG",2)
	W = VectorFunctionSpace(mesh,"CG",3)
	u = TrialFunction(V)
	v = TestFunction(V)
	class Boundary(SubDomain):
		def inside(self, x, on_boundary):
			return on_boundary
	u_e = Expression(("pi*x[0]*cos(pi*x[0]*x[1])","-pi*x[1]*cos(pi*x[0]*x[1])"))
	u_exact = project(u_e,W)
	f = Expression(("pi*pi*(pi*x[0]*(x[0]*x[0]+x[1]*x[1])*cos(pi*x[0]*x[1])+2*x[1]*sin(pi*x[0]*x[1]))",\
	                 " -pi*pi*(pi*x[1]*(x[0]*x[0]+x[1]*x[1])*cos(pi*x[0]*x[1])+2*x[0]*sin(pi*x[0]*x[1]))"))
	boundary = Boundary()
	bound = FacetFunction("size_t",mesh)
	boundary.mark(bound,1)
	#plot(bound,interactive=True)

	bc1 = DirichletBC(V, u_e, boundary)
	#bc2 = DirichletBC(Q, 0, boundary)
	bcs = [bc1]
	mu = Constant(mu)
	lamb = Constant(lamb)
	a = mu*inner(grad(u),grad(v))*dx - lamb*inner(grad(div(u)),v)*dx
	L = inner(f,v)*dx
	u_ = Function(V)
	solve(a==L,u_,bcs)
	return u_exact, u_, mesh




N = [8,16,32,64]
mu = 1.0
b = np.zeros(len(N))
a = np.zeros(len(N))
for lamb in[1,100,10000]:
	for i in range(len(N)):
		u ,u_exact,mesh = solve_me(N[i],mu,lamb)
		e = errornorm(u_exact, u, norm_type = "l2",degree_rise = 2)
		#print "N : ",N ,"lambda: ",lamb
		print "N: %e Lambda: %e Errornorm: %e"%(N[i],lamb ,e)
		#plot(u_exact)
		#plot(u)
		#interactive()
		b[i] = np.log(e)#/u_norm)
		a[i] = np.log(mesh.hmin())#1./(N[i]))
		#print u_norm
		#u_norm = sum(u_ex*u_ex)**0.5
	plt.plot(a,b, label =("lamb %.f "%(lamb)))
	axes = plt.gca()
	legend = axes.legend(loc='lower right', shadow=True)
	A = np.vstack([a,np.ones(len(a))]).T
	alpha,c = np.linalg.lstsq(A,b)[0]
	print "lambda: ",lamb
	print "alpha: ", alpha
	print "c: ", c
	print "-----------"
plt.show()
