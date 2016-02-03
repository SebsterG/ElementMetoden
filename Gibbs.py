import numpy as np 
import matplotlib.pyplot as plt


def func(x):
	return np.logical_and(x>= 1./3. , x<2./3.)
	
def u(N,K,x):
	b = np.zeros(K)
	dx = 1./N
	for i in range(K):
		b[i] = 2 *np.trapz(func(x)*np.sin(i*np.pi*x))*dx
	u = 0
	for j in range(K):
		u += b[j]*np.sin(j*np.pi*x)
	return u

def L_2(N,K):
	x = np.linspace(0,1,N+1)
	diff = u(N,K,x)-func(x)
	return np.sqrt(sum(diff**2))/len(diff)

def sup(N,K):
	x = np.linspace(0,1,N+1)
	return max(u(N,K,x) - func(x))

def sobolev(N,K):
	x = np.linspace(0,1,N+1)
	diff = u(N,K,x)-func(x)
	deriv = (diff[1:] - diff[:-1])/ (x[1]-x[0])
	return np.sqrt(sum(deriv**2.))/len(deriv) + np.sqrt(sum(diff**2.))/len(diff)
for i in range(1,100):
	i += 1000
	print sobolev(2*i,100)


#plt.plot(u(1000,100))f
#plt.show()


#plt.plot(b[1]*np.sin(2*np.pi*x_1)+b[2]*np.sin(3*np.pi*x_1))
#plt.show()


















