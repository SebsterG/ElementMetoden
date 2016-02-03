import numpy as np
import matplotlib.pyplot as plt
import sympy as sym


#def integral_symbols():
x, h = sym.symbols('x h')
#print sym.integrate(sym.cos(x),(x,0,sym.pi/2.0),(h,2))
up = (1./h)*x-(1./h)*(0.5-h)
down = (-1./h)*x+(1./h)*(0.5+h)
first = sym.integrate(up**2 ,(x,0.5-h,0.5)) 
last = sym.integrate(down**2,(x,0.5,0.5+h))
Sum = first.subs(h,0.1) + last.subs(h,0.1)

first_1 = sym.integrate(up**2+(sym.diff(up,x))**2 ,(x,0.5-h,0.5)) 
last_1 = sym.integrate(down**2+(sym.diff(down,x))**2 ,(x,0.5,0.5+h)) 
Sum_1 = first_1.subs(h,0.1) + last_1.subs(h,0.1) 


print "L_2 Analytical solution: ",sym.sqrt(Sum)
print "H_1 Analytical solution: ",sym.sqrt(Sum_1) 





def hat(x,h):
	if x>(0.5-h) and x<0.5:
		return (1./h)*x-(1./h)*(0.5-h)
	if x>=0.5 and x<(0.5+h):
		return -(1./h)*x+(1./h)*(0.5+h)
	else:
		return 0.0	
def H_1(N,u):
	dx = 1.0/(N-1)
	deriv = (u[1:] - u[:-1])/dx
	return np.sqrt(sum(deriv**2.)/len(deriv)) + np.sqrt(sum(u**2.)/len(u))

def L_2(u):
	return np.sqrt(sum(u**2.)/len(u))

def y(N,h):
	x = np.linspace(0,1.0,N+1)
	y = np.zeros(N+1)
	for i in range(len(x)):
		y[i] = hat(x[i],h)
	return y
h = 0.1
for N in (10,100,1000,10000):
	print"L_2 = %.4f , N = %.4f" %(L_2(y(N,h)),N)
	print"H_1 = %.4f , N = %.4f" %(H_1(N,y(N,h)),N)
	print"---------------------------------"






