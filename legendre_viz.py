###VISUALIZATION FOR 1st 10 FOURIER TRANSFORMS & FIRST 10 LEGENDRE POLYNOMIALS

#Package Imports
import numpy as np
import matplotlib.pylab as plt
from scipy.special import legendre

#Plotting first 10 Fourier expansions
#Define parameters
nmax = 50
nx = 50
x = np.linspace(0,2*np.pi,nx)
y = np.zeros((nmax,nx))

#Fourier transform
for n in range(1,nmax):
    xn = (-2)*np.sin(n*x)/n
    y[n,:] = y[n-1,:]+xn

#Loop Fourier transform 10 times
for i in range(1,11):
    plt.plot(x, y[i,:]+1, label= 'n =' + str(i))

#Plotting details
plt.title("Fourier expansions 1-10 for f(x) = x")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.show()

#Plotting first 10 Legendre Polynomials

#Define parameters
n = 1000
xv = np.linspace(-1,1,n)

#Define a Legenre func using scipy
def Legendre(x,N):
    leg = legendre(N)
    P_N = leg(x)
    return P_N

#Loop function for 10 iterations
for i in range(1,11):
    func = Legendre(xv, i)
    plt.plot(xv, func, label= 'n =' + str(i))

#Plotting details
plt.title('Legendre polynomials 1-10 for f(x) = x')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.show()
