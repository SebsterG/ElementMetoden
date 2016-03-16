from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
set_log_active(False)

def solver_1(N ,mu, beta):
    #print exp(1000)
    mu = Constant(mu)
    mesh = UnitSquareMesh(N, N)

    V = FunctionSpace(mesh, "CG", 1)
    V2 = FunctionSpace(mesh, "CG",2)

    class Right(SubDomain):
    	def inside(self, x, on_boundary):
    		return on_boundary and near(x[0],1)

    class Left(SubDomain):
    	def inside(self, x, on_boundary):
    		return on_boundary and near(x[0],0)


    right = Right()
    left = Left()
    bound = FacetFunction("size_t", mesh)
    bound.set_all(0)
    right.mark(bound, 1)
    left.mark(bound,2)
    #plot(bound); interactive()

    bc0 = DirichletBC(V, 1, right)
    bc1 = DirichletBC(V, 0, left)
    bcs = [bc0,bc1]

    u = TrialFunction(V)
    v = TestFunction(V)
    v = v + beta*mesh.hmin()*v.dx(0)
    f = Constant(0)

    u_exact = interpolate(Expression("(exp(x[0]/mu)-1)/(exp(1/mu)-1)",mu=mu),V2)
    # sette 0 overalt utenom paa 1
    #print "norm of exact: ", norm(u_exact,"l2")
    #plot(u_exact);interactive()

    a = mu*inner(grad(u), grad(v))*dx + u.dx(0)*v*dx
    F = f*v*dx

    u_1 = Function(V)
    solve(a==F,u_1,bcs)
    #print project(u_exact - u_1,V)
    #print "mu = %.e , N = %.f Numerical Error: %.4e " %(mu,N,errornorm(u_exact, u_1, norm_type="l2",degree_rise = 4))

    #plot(interpolate(u_exact,V))
    #plot(u_1); interactive()
    #plt.show()
    #plt.savefig("Galerkin.png")
    return u_1,u_exact,V,mesh
#solver_1(64,0.0001,beta = 1.0)



N = [2,4,8,16,32,64]
mu = 0.01
b = np.zeros(len(N))
a = np.zeros(len(N))


for i in range(len(N)):
    u ,u_exact,V,mesh = solver_1(N[i],mu,beta=1.0)
    #ex_norm = assemble((u_exact.dx(0) - u.dx(0))**2*dx)
    #ey_norm = assemble((u_exact.dx(1) - u.dx(1))**2*dx)

    #e = sqrt(mesh.hmin()*ex_norm + mu*(ex_norm+ey_norm))
    e = errornorm(u_exact, u, norm_type = "l2",degree_rise = 3)
    print e
    #u_norm = norm(u_exact,"H1")
    #e = errornorm(u_exact,u,norm_type = "l2",degree_rise = 3)
    b[i] = np.log(e)#/u_norm)
    a[i] = np.log(mesh.hmin())#1./(N[i]))
    #print u_norm
    #u_norm = sum(u_ex*u_ex)**0.5
A = np.vstack([a,np.ones(len(a))]).T
alpha_1,c_1 = np.linalg.lstsq(A,b)[0]

print "mu = ",mu, "C =  ", exp(c_1), "alpha =  ", alpha_1

plt.plot(a,b,'o',label='Original data',markersize = 10)
plt.plot(a,alpha_1*a+c_1,'r',label='Fitted line')
plt.legend();plt.show()

"""

mu =[1, 0.1 ,0.01, 0.001]; N = [8,16,32,64];
for m in mu:
    for i in N:
        solver_1(i,m,beta=0.15)"""
