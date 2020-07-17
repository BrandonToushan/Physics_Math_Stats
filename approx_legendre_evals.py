####CALCULATING LEGENDRE POLYNOMIALS via MATRIX CALCULUS & FINITE DIFF APPROXMIATION 
####& ORTHONORMALIZING EIGENVECTORS via GRAM-SCHMIDT 

#WIKIPEDIA: "Legendre polynomials (named after Adrien-Marie Legendre, who discovered them in 1782) are a system of complete and orthogonal polynomials, 
#with a vast number of mathematical properties, and numerous applications. They can be defined in many ways, and the various definitions highlight 
#different aspects as well as suggest generalizations and connections to different mathematical structures and physical and numerical applications.

#Package Imports
import numpy as np
import math as m
import matplotlib.pylab as plt
%matplotlib inline

#Legendre Equation done in Python:

#Construct a matrix that approximates the deriv

#Define range and fineness of function
xmin = (-1)
xmax = (1)
nn = 300
dx = (xmax-xmin)/nn
x = np.linspace(xmin,xmax,nn,True)

#Define $d/dx$
d_x = np.zeros([nn,nn],dtype=float)
d_x[0,0]=-1.0/dx
d_x[0,1]=1.0/dx
d_x[-1,-1]=1.0/dx
d_x[-1,-2]=-1.0/dx

for i in range(1,nn-1) :
    d_x[i,i-1]=-0.5/dx
    d_x[i,i+1]=0.5/dx
    
#Define the function we want to differentiate
fx = np.array([x_i for x_i in x])

#Finally, we compute the derivative
dfx_dx = d_x.dot(fx)

#Verify code still correctly computes deriv of f(x) = x

#Plot, f(x) is red, f'(x) is blue
plt.figure()
plt.title("f(x), df(x)/dx vs x")
plt.ylabel('f(x),df(x)/dx')
plt.xlabel('x')
plt.plot(x,fx,'rs')
plt.plot(x,dfx_dx,'gs')

#Repeat for Various polynomials

#i) f(x) = x^2
fx = [x_i**2 for x_i in x]
dfx_dx = d_x.dot(fx)

#plt.figure()
#plt.ylabel('f(x),df(x)/dx')
#plt.xlabel('x')
#plt.plot(x,fx,'rs')
#plt.plot(x,dfx_dx,'bo')

#ii) f(x) = exp(x)
fx = [np.exp(x_i) for x_i in x]
dfx_dx = d_x.dot(fx)

#plt.figure()
#plt.ylabel('f(x),df(x)/dx')
#plt.xlabel('x')
#plt.plot(x,fx,'rs')
#plt.plot(x,dfx_dx,'bo')

#iii) f(x) = sin(2pix)
fx = np.array([np.sin(2*np.pi*x_i) for x_i in x])
dfx_dx = d_x.dot(fx)

#plt.figure()
#plt.ylabel('f(x),df(x)/dx')
#plt.xlabel('x')
#plt.plot(x,fx,'rs')
#plt.plot(x,dfx_dx,'bo')

# Write down the algebraic expression for a finite diff approx:

#Define $d/dx$
d_x = np.zeros([nn,nn],dtype=float)
d_x[0,0]=-1.0/dx
d_x[0,1]=1.0/dx
d_x[-1,-1]=1.0/dx
d_x[-1,-2]=-1.0/dx

for i in range(1,nn-1) :
    d_x[i,i-1]=-0.5/dx
    d_x[i,i+1]=0.5/dx

#Take derivative of d_x to get d2_x 
d2_x = d_x.dot(d_x)

#Check that d2_x takes the 2nd derivative by plotting

#try f(x) = x^2
fx = [x_i**2 for x_i in x]
d2fx_dx = d2_x.dot(fx)
plt.figure()
plt.ylabel('f(x),d^2f(x)/dx^2')
plt.xlabel('x')
plt.xlim(-.98,.98)
plt.plot(x,fx, linewidth = 5, color = 'red')
plt.plot(x,d2fx_dx, linewidth = 5, color = 'green')

#try f(x) = x^3
fx = [x_i**3 for x_i in x]
d2fx_dx = d2_x.dot(fx)
plt.figure()
plt.ylabel('f(x),d^2f(x)/dx^2')
plt.xlabel('x')
plt.plot(x,fx, linewidth = 5, color = 'red')
plt.plot(x,d2fx_dx,linewidth =5, color = 'green' )
plt.xlim(-.98,.98)

#Therefore d2_x works as a second derivative function

#Construct matrices to represent multiplication by functions

#i) (1-x^2)
C1 = np.zeros([nn,nn],dtype = float)
for i in range(0,nn): 
    C1[i,i]  = 1 - x[i]**2

#ii) (2x)
C2 = np.zeros([nn,nn], dtype = float)

for i in range(0,nn):
    C2[i,i] = -2*x[i]

#iii) Plot C1 and C2 x f(x) = 1 to check
fx = np.array([1 for x_i in x])
fx = np.zeros([nn,1],dtype = float)
for i in range(0,nn):
    fx[i,0] = 1

C1_fx = C1.dot(fx)
C2_fx = C2.dot(fx)

#Check functions by plotting 

#i) (1-x^2)
plt.figure()
plt.title('(1-x^2) vs x')
plt.ylabel('(1-x^2)')
plt.xlabel('x')
plt.plot(x,C1_fx, linewidth = 5, color= 'green')

#ii) (-2x)
plt.figure()
plt.title('-2x vs x')
plt.ylabel('(-2x)')
plt.xlabel('x')
plt.plot(x,C2_fx, linewidth = 5, color = 'red')

Define a matrix L that approximates the Legendre Equation
L =  C1.dot(d2_x) + C2.dot(d_x)

#Check 
print(L)
print(L.shape)

#L is a matrix, made up of d_x,d2_x,C1 and C2

#L acting on f(x) = x
fx = np.array([x_i for x_i in x])
L_fx = L.dot(fx)

#transpose of L acting on f(x) = x
fx = np.array([x_i for x_i in x])
L_trans = L.transpose()
L_trans_fx = L_trans.dot(fx)

#plotting
plt.figure()

plt.title('L_fx,L_trans_fx vs x')
plt.ylabel('L_fx')
plt.xlabel('x')
plt.xlim(-0.95,0.95)
plt.ylim(-2,2)
plt.plot(x,L_fx, color= 'green')
plt.plot(x,L_trans_fx, color= 'red')

Compare L_fx and L_trans_fx for various functs

#i) f(x) = x^3
fx = np.array([x_i**3 for x_i in x])
L_fx = L.dot(fx)

#transpose of L acting on f(x) = x^3
L_trans = L.transpose()
L_trans_fx = L_trans.dot(fx)

#plotting
#plt.figure()
#plt.xlim(-.95,.95)
#plt.ylim(-3,3)
#plt.title('L_fx,L_trans_fx vs x')
#plt.ylabel('L_fx')
#plt.xlabel('x')
#plt.plot(x,L_fx, color= 'purple')
#plt.plot(x,L_trans_fx, color= 'yellow')

#ii) f(x) = exp(x)
fx = np.array([np.exp(x_i) for x_i in x])
L_fx = L.dot(fx)

#transpose of L acting on f(x) = exp(x)
L_trans = L.transpose()
L_trans_fx = L_trans.dot(fx)

#plotting
#plt.figure()
#plt.xlim(-.95,.95)
#plt.title('L_fx,L_trans_fx vs x')
#plt.ylabel('L_fx')
#plt.xlabel('x')
#plt.ylim(-3,3)
#plt.plot(x,L_fx, color= 'purple')
#plt.plot(x,L_trans_fx, color= 'yellow')

#iii) f(x) = sin(x)
fx = np.array([np.sin(x_i) for x_i in x])
L_fx = L.dot(fx)

#transpose of L acting on f(x) = sin(x)
L_trans = L.transpose()
L_trans_fx = L_trans.dot(fx)

#plotting
#plt.figure()
#plt.xlim(-.95,.95)
#plt.ylim(-3,3)
#plt.title('L_fx,L_trans_fx vs x')
#plt.ylabel('L_fx')
#plt.xlabel('x')
#plt.plot(x,L_fx, color= 'purple')
#plt.plot(x,L_trans_fx, color= 'yellow')

#Calculating the eigenvalues of L

#import numpy's lin alg package to use for calculations
import numpy.linalg 
import scipy.linalg

#i) Calculate eigenvalues of real, symmetric matrix L
e_vals, e_vecs = scipy.linalg.eigh(L)

#ii) Arranging eigenvalues in decreasing order
sorted_evals = sorted(e_vals, reverse = True)
n = np.arange(0,300,1)
n_n = (-1)*n
n_2n = (n+1)*n_n

#plotting
plt.figure()
plt.title('Lambda_n vs n and -n(n+1) vs n')
plt.ylabel('Lambda_n')
plt.xlabel('n')
plt.scatter(n,sorted_evals, color = 'red')
plt.scatter(n,n_2n, color = 'green')

#Finding the eigenvectors of L

#From above:
e_vals, e_vecs = scipy.linalg.eig(L)

#Find poistion of (almost) 0 value in list so as to find corresponding eigenvec

#Convert np array to list and then round in order to use index func to find position of "0" eigenvec"
e_val_list = e_vals.tolist()

#As index(0) doesn't work due to the value being not quite 0, the value had to be found using brute force
#print(e_val_list)
a = e_val_list.index(-7.213923634679477e-13+0j)

#Thus the eigenvec corresponding to the 0 eigenval is
P_0 = e_vecs[a]
#Checking eigenvec
L_P_0 = L.dot(P_0)

#Plotting
plt.figure()
plt.title('L_P_0 vs x ')
plt.ylabel('L_P_0')
plt.xlabel('x')
plt.plot(x,L_P_0, color= 'green')

#Other eigenvectors via Gram-Schmidt

#i) Trial function (DNS):

t_1  = np.array([x_i for x_i in x])

##ii) Putative eigenvector (DNS):|P_1>, via Gram-Schmidt

#Gram-Schmidt projecting |P_0> orthogonally onto the line spanned by t
a = np.inner(P_0,t_1)
b = np.inner(t_1,t_1)
P_1 = t_1*(a*(b))

#iii) Plotting (Submit):

#L acting on P_1
L_P1 = L.dot(P_1)

plt.figure()
plt.title('t_1,P_1 and L_P_1 vs x ')
plt.xlabel('x')
plt.plot(x,t_1, color = 'red')
plt.plot(x,P_1, color = 'blue')
plt.plot(x,L_P1, color= 'g')

#Is P_1 an eigenvector of L?

#For P_1 to be a an eigenvector of L  L|P_1> = lambda|P_1>

#Ideally should find a matrix of one constant (the eigenval)

e_val_1 = P_1/L

#Thus the corresponding e_val is 0.03006337

e_val_1 = 0.03006337

#i) x^2
t_2  = np.array([x_i**2 for x_i in x])

#Gram-Schmidt projecting |P_1> orthogonally onto the line spanned by t_2
a = np.inner(P_0,t_2)
b = np.inner(t_2,t_2)
P_2 = t_2*(a*(b))

#ii) x^3
t_3 = np.array([x_i**3 for x_i in x])

#Gram-Schmidt projecting |P_2> orthogonally onto the line spanned by t_3
a = np.inner(P_0,t_3)
b = np.inner(t_3,t_3)
P_3 = t_3*(a*(b))

#iii) x^3
t_4 = np.array([x_i**4 for x_i in x])

#Gram-Schmidt projecting |P_2> orthogonally onto the line spanned by t_3
a = np.inner(P_0,t_3)
b = np.inner(t_3,t_3)
P_4 = t_3*(a*(b))

#Plotting

#i) |P_2>
L_P_2 = L.dot(P_2)
plt.figure()
plt.title('t_2,P_2 and L_P_2 vs x ')
plt.xlabel('x')
plt.plot(x,t_2, color = 'r')
plt.plot(x,P_2, color = 'b')
plt.plot(x,L_P_2, color= 'g')

#ii) |P_3>
L_P_3 = L.dot(P_3)
plt.figure()
plt.title('t_3,P_3 and L_P_3 vs x ')
plt.xlabel('x')
plt.plot(x,t_3, color = 'r')
plt.plot(x,P_3, color = 'b')
plt.plot(x,L_P_3, color= 'g')

#ii) |P_4>
L_P_4 = L.dot(P_4)
plt.figure()
plt.title('t_4,P_4 and L_P_4 vs x ')
plt.xlabel('x')
plt.plot(x,t_4, color = 'r')
plt.plot(x,P_4, color = 'b')
plt.plot(x,L_P_4, color= 'g')

#Finding Eigenvalues

#i) |P_2>

#For P_2 to be a an eigenvector of L  L|P_2> = lambda|P_2>
#Ideally should find a matrix of one constant (the eigenval)

e_val_2 = P_2/L
#print(e_val_2)

#Thus the corresponding e_val is 0.16272140
e_val_2 = 0.16272140

#ii) |P_3>

#For P_3 to be a an eigenvector of L  L|P_3> = lambda|P_3>
#Ideally should find a matrix of one constant (the eigenval)

e_val_3 = P_3/L
#print(e_val_3)

#Thus the corresponding e_val is 0.01425235
e_val_3 = 0.01425235

#iii) |P_4>
#For P_4 to be a an eigenvector of L  L|P_4> = lambda|P_3>
#Ideally should find a matrix of one constant (the eigenval)

e_val_4 = P_4/L
#print(e_val_4)

#Thus the corresponding e_val is 0.00717384
e_val_4 = 0.00717384
