from dolfin import *
import time
set_log_active(False)

mesh = Mesh("cylpipe.xml")
#plot(mesh,interactive=True)

V = VectorFunctionSpace(mesh,"CG",2)
Q = FunctionSpace(mesh,"CG",1)

u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

class Outlet(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[2],1.0)
class Inlet(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary and near(x[2],0.0)
class Nos(SubDomain):
	def inside(self,x,on_boundary):
		return on_boundary

outlet = Outlet()
inlet = Inlet()
nos = Nos()
bound = FacetFunction("size_t", mesh)
bound.set_all(0)
nos.mark(bound, 1)
outlet.mark(bound, 2)
inlet.mark(bound, 3)
#plot(bound, interactive = True)
n = FacetNormal(mesh)
a = 1.0
b = 1.0
#p_0 = Expression(" 1 + a*sin(pi*t) + b*sin(4*pi*t)",t=0.0,a=a,b=b)
p_0 = Constant(1.0)
p_1 = Constant(0.0)

bc1 = DirichletBC(Q, p_0 , inlet)
bc2 = DirichletBC(V, (0.0,0.0,0.0), nos)
bc3 = DirichletBC(Q, p_1, outlet)
bcs = [bc2]
bcp = [bc1,bc3]
u0 = Function(V)
u1 = Function(V)
u_star = Function(V)
p0 = Function(Q)
p1 = Function(Q)

dt = 0.0001
rho = Constant(1000)
mu = Constant(0.001002)
K = Constant(dt)
def sigma(u, p):
    return 2.0*mu*sym(grad(u))-p*Identity(3)#(len(u))
def eps(u):
    return sym(grad(u))

f = Constant((0.0,0.0,0.0))

F1= (1.0/K)*rho*inner(u-u0,v)*dx + \
inner(grad(u0)*u0, v)*dx+\
inner(sigma(0.5*(u+u0),p0),eps(v))*dx \
- mu*dot(dot(grad(0.5*(u+u0)),n),v)*ds(3) + \
dot(p0*n,v)*ds(3) - inner(f,v)*dx

a2 = K *inner(grad(p),grad(q))*dx
L2 = K *inner(grad(p0),grad(q))*dx - rho*inner(div(u_star),q)*dx

a3 = rho*inner(u,v)*dx
L3 = rho*inner(u_star,v)*dx - K*inner(grad(p1-p0),v)*dx
T = 0.1
t = dt
counter = 0
import shutil
shutil.rmtree('results')
file_u = File("results/velocity.pvd")
file_p = File("results/pressure.pvd")
while t < T + DOLFIN_EPS:
    start_time = time.clock()
    # Update pressure boundary condition
    #v_theta.t = t
    p_0.t = t
    solve(lhs(F1)==rhs(F1),u_star,bcs)

    #pressure correction
    solve(a2==L2,p1,bcp)
    #print norm(p1)

    #last step
    solve(a3==L3,u1)#,bcs)

    u0.assign(u1)
    p0.assign(p1)
    t += dt
    print "norm",norm(u1,"l2")
    print "Timestep: ", t
    file_u << u1
    file_p << p1
    print time.clock() - start_time, "seconds"

    #plot(u1,rescale=True)
#print assemble(u1[0]*dx)/assemble(1.0*dx(mesh))
#interactive()
