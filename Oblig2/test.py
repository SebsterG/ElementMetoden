from sympy import *
x, y, pi = symbols('x y pi')
#print ((x+y)**2 * (x+1)).expand()
#(pi*x*cos(pi*x*y), -pi*y*cos(pi*x*y))
f1 = simplify(diff(pi*x*cos(pi*x*y), x,x) + diff(pi*x*cos(pi*x*y), y,y))

f2 = simplify(diff(-pi*y*cos(pi*x*y), x,x) + diff(-pi*y*cos(pi*x*y), y,y))

f = lambdify((x, y), [f1, f2])
print f

#k, m, n = symbols('k m n', integer=True)
#f, g, h = map(Function, 'fgh')
