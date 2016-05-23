from dolfin import *
import time
import matplotlib.pyplot as plt
jacobi_time =[]
"""parameters["krylov_solver"]["relative_tolerance"] = 1.0e-8
parameters["krylov_solver"]["absolute_tolerance"] = 1.0e-8
parameters["krylov_solver"]["monitor_convergence"] = False
parameters["krylov_solver"]["report"] = False
parameters["krylov_solver"]["maximum_iterations"] = 50"""
def solve_time(A1,b1,bcs,pre,type1):
    t0 = time.time()
    u1 = Function(V)

    [bc.apply(A1,b1) for bc in bcs]
    pc = PETScPreconditioner(pre)
    sol = PETScKrylovSolver(type1,pc)
    prm = sol.parameters
    prm["maximum_iterations"] = 100000
    prm["relative_tolerance"] = 1.0e-10
    prm["absolute_tolerance"] = 1.0e-10
    prm["monitor_convergence"] = True
    #sol.set_tolerances(1.0e-8, 1.0e-8,10000)
    sol.solve(A1,u1.vector(),b1)
    #[bc.apply(A1,b1) for bc in bcs]
    #solve(A1, u1, b1)#, "cg", "hypre_amg")
    #solve(A1 == b1, u1, bcs=bcs)#,solver_parameters={"linear_solver": "default"},form_compiler_parameters={"optimize": True})
    t1= time.time()
    return t1-t0
for N in [32,64,128,256,512,1024]:
    mesh = UnitIntervalMesh(N)

    V = FunctionSpace(mesh,"CG",1)
    u = TrialFunction(V)
    v = TestFunction(V)

    class Boundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary
    boundary = Boundary()
    bound = FacetFunction("size_t",mesh)
    boundary.mark(bound,1)
    bc0 = DirichletBC(V,0,boundary)
    bcs = [bc0]

    #f = Expression("x[0]")
    f = matrix(random.random(N)).transpose()

    a = inner(grad(u),grad(v))*dx
    L = f*v*dx

    A1 = assemble(a)
    b1 = assemble(L)
    #print solve_time(a,L,bcs)

    t1 = solve_time(A1,b1,bcs,"jacobi","cg")
    print t1
    jacobi_time.append(t1)
import pylab
pylab.plot(jacobi_time)
pylab.show()
