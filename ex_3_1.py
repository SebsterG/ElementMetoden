import numpy as np 


def u_random(N):
	return np.random.random_sample((N,))
def u_sin(N,k):
	x = np.linspace(0,1,N+1)
	return np.sin(k*np.pi*x)

def H_1(N,u):
	dx = 1.0/(N-1)
	deriv = (u[1:] - u[:-1])/ dx
	return np.sqrt(sum(deriv**2.)/len(deriv)) + np.sqrt(sum(u**2.)/len(u))

def L_2(N,u):
	return np.sqrt(sum(u**2.)/len(u))

for N in (10,100,1000,10000):
	print "L_2: ",L_2(N,u_random(N)), " N = ",N
	print "H_1: ",H_1(N,u_random(N)), " N = ",N
for k in (1,10,100):
	print "L_2 = %.2f k = %.2f, N = %.2f " %(L_2(N,u_sin(100,k)),k,N)
	print "H_1 = %.2f k = %.2f, N = %.2f " %(H_1(N,u_sin(100,k)),k,N)