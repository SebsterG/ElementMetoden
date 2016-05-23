from dolfin import *
from numpy import *
import time
lu_time =[]; ilu_time = []; jacobi_time=[]; amg_time=[];
Ns = [32,64,128,256,512,1024]
parameters["krylov_solver"]["relative_tolerance"] = 1.0e-8
parameters["krylov_solver"]["absolute_tolerance"] = 1.0e-8
parameters["krylov_solver"]["monitor_convergence"] = False
parameters["krylov_solver"]["report"] = False
parameters["krylov_solver"]["maximum_iterations"] = 50000

def solving_time(A1,b1,solver):
    U = Function(V)
    t0 = time.time()
    bcs =[]
    if len(solver) == 2:
        solve(A,U.vector(),b,solver[0],solver[1]);
        """[bc.apply(A1,b1) for bc in bcs]
        pc = PETScPreconditioner(solver[1])
        sol = PETScKrylovSolver(solver[0],pc)
        prm = sol.parameters
        prm["maximum_iterations"] = 100000
        prm["relative_tolerance"] = 1.0e-10
        prm["absolute_tolerance"] = 1.0e-10
        prm["monitor_convergence"] = True
        #sol.set_tolerances(1.0e-8, 1.0e-8,10000)
        sol.solve(A1,U.vector(),b1)"""
    else:
        solve(A,U.vector(),b,solver[0]);
        """[bc.apply(A1,b1) for bc in bcs]
        pc = PETScPreconditioner("default")
        sol = PETScKrylovSolver(solver[0],pc)
        prm = sol.parameters
        prm["maximum_iterations"] = 100000
        prm["relative_tolerance"] = 1.0e-10
        prm["absolute_tolerance"] = 1.0e-10
        prm["monitor_convergence"] = True
        sol.solve(A1,U.vector(),b1)"""

    t1 = time.time()
    return t1-t0
for N in Ns:
    mesh = UnitSquareMesh(N,N)
    print "N", N, "Dofs: ", mesh.num_vertices()
    V = FunctionSpace(mesh,"Lagrange",1)
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Expression("sin(x[0])")
    #f = matrix(random.random(N)).transpose()
    a = inner(grad(u),grad(v))*dx
    L = f*v*dx
    A = assemble(a)
    b = assemble(L)
    t2 = solving_time(A,b,["lu"])
    print"Time for lu", t2
    lu_time.append(t2)
    #t2= solving_time(A, b,["cg","jacobi"])
    #print"Time for cg", t2
    #jacobi_time.append(t2)
    """t2 = solving_time(A, b,["cg","ilu"])
    print"Time for cg/ilu", t2
    ilu_time.append(t2)
    t2 = solving_time(A, b,["cg","amg"])
    print"Time for cg/amg",t2
    amg_time.append(t2)"""

import pylab
pylab.plot(Ns,lu_time)
pylab.plot(Ns,ilu_time)
pylab.plot(Ns,jacobi_time)
pylab.plot(Ns,amg_time)
pylab.xlabel("unknowns")
pylab.ylabel("Time/sec")
pylab.legens(["lu","ilu","jacobi","amg"])
pylab.show()
