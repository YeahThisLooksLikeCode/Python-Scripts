import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

########### boundary condition functions ########
def boundaries(x,xmin,xmax,Nx): 
	dx = (xmax - xmin)/Nx
	for i in range(Nx - 1):
		x[i] = xmin + (i + 0.5)*dx
	return x

def f(x):    #u(x,0) = f(x)
    return np.sin(np.pi * x)

def g(t):    #u(minx,t) = g(t)
    return f(0.)
def h(t):    #u(maxx,t) = h(t)
    return f(1.)

def u(x): #u(x,0) set to zero for first half and one for second half.
	for i in range(Nx-1):
		if x[i] <= 0.5:
			return 0.
		else:
			return 1.
			
def w(v,Nx):
	w = 0.
	for i in range(Nx - 1):
		w += v[i,:]
	return w
	
########### Functions ############
def hybridSolver(mu,v):
    a,b,c, d = getAbcTheta(Nx)
    theta = 0.5
    for n in range(Nt - 1):

        d[0] = (1. - (1. - theta)*mu)*v[0,n] + (1.-theta)*mu*v[1,n] # set first step, i =0
        d[1:-1] = (1. - 2.*(1. - theta)*mu)*v[1:-1,n] + (1.-theta)*mu * (v[:-2,n] + v[2:,n]) #otherwise
        d[-1] = (1.-theta)* mu * v[-2,n] + (1. - (1. - theta)*mu)*v[-1,n] # set last step, i = Nx -1
        v[:,n+1] = thomasAlgorithm(a,b,c,d) 
    return v

def thomasAlgorithm(a,b,c,d):
    N = len(d)
    x = np.empty(N)
    y = np.empty(N)
    y[0] = c[0]/b[0]
    x[0] = d[0]/b[0] 
    for j in range(1, N):
        y[j] = c[j] / (b[j] - a[j]*y[j-1])
        x[j] = (d[j] + a[j]*x[j-1]) / (b[j] - a[j]*y[j-1])
    for j in reversed(range(N-1)):
        x[j] +=  y[j]*x[j+1]

    return x

def getAbcTheta(Nx): 
    theta = 0.5
    n= (Nx - 2)
    a = np.ones(Nx) * theta * mu # set a as theta*mu at all points
    b = np.ones(Nx)*(1.+theta*mu) # b equals 1 + theta*mu at first step and last
    b[1:-1] += 2*theta * mu # and 1+ 2*theta*mu elsewhere
    c = a.copy() # set c as theta*mu at all points execpt the last Nx - 1
    c[-1] = 0  # i = Nx - 1 = 0
    a[0] = 0 # and rewrite the first value of a as zero
    d = np.empty(Nx) 
    return a,b,c,d

##### Set variables and run functions #####
mu = 0.5
xmin= 0.
xmax = 1.0
Nx = 20
Nt = 50
dt = 0.01
t = np.arange(Nt)*dt
x = np.linspace(xmin, xmax, Nx)

#### for boundary conditions set to a sine function #######
v = np.empty([Nx,Nt])
v[:,0] = f(x)
v[0,:] = g(t)
v[-1, :] = h(t)

v = hybridSolver(mu,v)
print (v)

fig = plt.figure () 
ax = fig.add_subplot(111 , projection ='3d')
T, X = np.meshgrid (t, x)
ax.plot_wireframe (X, T, v)

ax.set_xlabel(r'$x_i$')
ax.set_ylabel(r'$t^{<n>}$')
ax.set_zlabel(r'$v i^{<n>}$')
plt.show ()

######### for boundary conditions set to function u(x,0)= either 0 or 1 ######## 
v = np.empty([Nx,Nt])
v[:,0] = u(x)
v[0,:] = g(t)
v[-1, :] = h(t)

v = hybridSolver(mu,v)
print (v)

fig = plt.figure ()
ax = fig.add_subplot(111 , projection ='3d')
T, X = np.meshgrid (t, x)
ax.plot_wireframe (X, T, v)

ax.set_xlabel(r'$x_i$')
ax.set_ylabel(r'$t^{<n>}$')
ax.set_zlabel(r'$v i^{<n>}$')
plt.show ()

w = w(v,Nx)
print (w)
plt.plot(w)
plt.show()
