import numpy as np 
import matplotlib.pyplot as plt

def hat(x):
	if x>=0.2 and x<0.3:
		return 10.0*x-2.0
	else:
		return -10.0*x + 4.0

def H_1(N,u):
	dx = 1.0/(N-1)
	deriv = (u[1:] - u[:-1])/dx
	return np.sqrt(sum(deriv**2.)/len(deriv)) + np.sqrt(sum(u**2.)/len(u))

def L_2(u):
	return np.sqrt(sum(u**2.)/len(u))
def y(N):
	x = np.linspace(0.2,0.4,N+1)
	y = np.zeros(N+1)
	for i in range(len(x)):
		y[i] = hat(x[i])
	#plt.plot(x,y)
	#plt.show()	
	return y

for N in (10,100,1000,10000):
	print"L_2 = %.2f , N = %.2f" %(L_2(y(N)),N)
	print"H_1 = %.2f , N = %.2f" %(H_1(N,y(N)),N)








