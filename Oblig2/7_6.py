from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
set_log_active(False)
def solve_this_shit(N,degree_u,degree_p):
    mesh = UnitSquareMesh(N,N)

    V = VectorFunctionSpace(mesh,"CG",degree_u)
    Q = FunctionSpace(mesh,"CG",degree_p)
    VQ = V*Q

    u,p = TrialFunctions(VQ)
    v,q = TestFunctions(VQ)


    class Up(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1],1)
    class Down(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1],0)
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0],0)
    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[0],1)



    u_e = Expression(("sin(pi*x[1])","cos(pi*x[0])"))
    p_e = Expression("sin(2*pi*x[0])")

    up = Up()
    down = Down()
    left = Left()
    right = Right()
    bound = FacetFunction("size_t",mesh)
    bound.set_all(0)
    up.mark(bound,0)
    down.mark(bound,1)
    left.mark(bound,2)
    right.mark(bound,3)
    #plot(bound,interactive=True)
    ds = Measure("ds",subdomain_data=bound)
    n = FacetNormal(mesh)
    ds = ds[bound]


    bc1 = DirichletBC(VQ.sub(0), u_e, up)
    bc2 = DirichletBC(VQ.sub(1), p_e , up)
    bc3 = DirichletBC(VQ.sub(0), u_e, down)
    bc4 = DirichletBC(VQ.sub(1), p_e , down)
    bc5 = DirichletBC(VQ.sub(0), u_e, left)
    bc6 = DirichletBC(VQ.sub(1), p_e , left)
    bc7 = DirichletBC(VQ.sub(0), u_e, right)
    bc8 = DirichletBC(VQ.sub(1), p_e , right)
    bcs = [bc1,bc2,bc3,bc4,bc5,bc6,bc7,bc8]

    f = Expression(("pi*pi*sin(pi*x[1])-2*pi*cos(2*pi*x[0])","pi*pi*cos(pi*x[0])"))
    a = inner(grad(u), grad(v))*dx + div(u)*q*dx + div(v)*p*dx - inner(f, v)*dx

    up = Function(VQ)
    solve(lhs(a)==rhs(a),up,bcs)
    u_,p_ = up.split()


    """
    #print
    plot(u_exact,interactive= True)
    plot(u_, interactive = True)
    plot(p_exact,interactive = True)
    plot(p_, interactive = True)"""


    R = VectorFunctionSpace(mesh, 'R', 0)
    c = TestFunction(R)
    tau = 0.5*(grad(u_)+grad(u_).T)
    n = FacetNormal(mesh)
    forces = assemble(dot(dot(tau, n), c)*ds(1)).array()
    print "x-direction = {}, y-direction = {}".format(*forces)
    #print forces
    #force = forces[0]
    """u_ex = project(u_e,V)
    R1 = VectorFunctionSpace(mesh, 'R', 0)
    c1 = TestFunction(R1)
    tau = 0.5*(grad(u_ex)+grad(u_ex).T)
    n = FacetNormal(mesh)
    forces_ex = assemble(dot(dot(tau, n), c1)*ds(1)).array()
    print "x-direction = {}, y-direction = {}".format(*forces_ex)"""

    #print forces_ex
    #force_exact = forces_ex[1]

    forces_ex=-0.5*(-2+pi)
    return u_,u_e,p_,p_e,mesh,V,Q#, forces[0],forces_ex

N = [4,8,16,32,64]
def convergence_force(N,degree_u,degree_p):
    b = np.zeros(len(N))
    a = np.zeros(len(N))
    d = np.zeros(len(N))

    for i in range(len(N)):
        u_,ue,p_,pe,mesh,V,Q,forces,forces_ex = solve_this_shit(N[i],degree_u,degree_p)
        e_u = abs(forces_ex-forces)
        #e_u = errornorm(forces_ex,forces, norm_type = "l2",degree_rise = 3)

        print "N:",N[i] ," force_error ",e_u
        print "------------------------------"
        a[i] = np.log(mesh.hmin())
        b[i] = np.log(e_u)

    plt.plot(a,b,label = "P%.f-P%.f"%(degree_u,degree_p))

    A = np.vstack([a,np.ones(len(a))]).T
    alpha,c = np.linalg.lstsq(A,b)[0]
    print "force-alpha: ", alpha
    print "force-c: ", c
    plt.show()
#convergence_force(N,4,3)


def convergence_e(N,degree_u,degree_p):
    b = np.zeros(len(N))
    a = np.zeros(len(N))
    d = np.zeros(len(N))

    for i in range(len(N)):
        u_,ue,p_,pe,mesh,V,Q = solve_this_shit(N[i],degree_u,degree_p)
        e_u = errornorm(ue, u_, norm_type = "h1",degree_rise = 3)
        e_p = errornorm(pe, p_, norm_type = "l2",degree_rise = 3)

        #print "N : ",N ,"lambda: ",lamb
        print "N:",N[i] ," u-errornorm: ",e_u
        print "N:",N[i] ," p-errornorm: ",e_p
        print "------------------------------"
        a[i] = np.log(mesh.hmin())
        b[i] = np.log(e_u)
        d[i] = np.log(e_p)

    plt.plot(a,b)
    #plt.plot(d,b)
    A = np.vstack([a,np.ones(len(a))]).T
    alpha,c = np.linalg.lstsq(A,b)[0]
    alpha_p,c_p = np.linalg.lstsq(A,d)[0]
    print "u-alpha: ", alpha
    print "u-c: ", c
    print "-----------"
    print "p-alpha: ", alpha_p
    print "p-c: ", c_p
    plt.show()
convergence_e(N,4,3)
