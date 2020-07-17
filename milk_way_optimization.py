###DETERMING MASS OF MILKY WAY GALAXY VIA POTENTIAL ENERGY FITTING & OPTIMIZATION
###GIVEN ROTATIONAL SPEED DATA FOR MILKY WAY, FIT THREE COMPONENTS WITH THE SAME FUNCTIONAL FORM OF THE POTENTIAL

#PACKAGE IMPORTS
import numpy as np
import math as m
import scipy
import matplotlib.pyplot as plt
from scipy import optimize


#EMPTY LISTS TO FILL WITH DATA
delta_V_kms = []
V_kms = []
r_kpc = []
delta_r_kpc = []

n = 2 #number of lines to ignore at the start of the file
w = 58 #last line to keep
i = 0 #tracker

#OPEN DATAFILE, APPEND COLUMNS TO LISTS ONE BY ONE
with open("/Users/BrandonToushan/Documents/PythonStuff/LogRC_data.dat", 'r') as file:   
    for line in file:
        if i <= n:
            i += 1  
        if i > n and i < w:
            i += 1
            row = line.split()
            delta_V_kms.append(row[-1:])
            r_kpc.append(row[0:1])
            delta_r_kpc.append(row[1:2])
            V_kms.append(row[2:3])    
    file.close()
    
#Turning lists to numpy arrays
delta_V_kms = np.array(delta_V_kms)
V_kms = np.array(V_kms)
r_kpc = np.array(r_kpc)
delta_r_kpc = np.array(delta_r_kpc)

#Split the values at E and calculate the value
delta_V_kms_1 = np.apply_along_axis(lambda a: (a[0].split('E')[0]),1,delta_V_kms)
delta_V_kms_2 = np.apply_along_axis(lambda a: (a[0].split('E')[1]),1,delta_V_kms)

V_kms_1 = np.apply_along_axis(lambda a: (a[0].split('E')[0]),1,V_kms)
V_kms_2 = np.apply_along_axis(lambda a: (a[0].split('E')[1]),1,V_kms)

delta_r_kpc_1 = np.apply_along_axis(lambda a: (a[0].split('E')[0]),1,delta_r_kpc)
delta_r_kpc_2 = np.apply_along_axis(lambda a: (a[0].split('E')[1]),1,delta_r_kpc)

r_kpc_1 = np.apply_along_axis(lambda a: (a[0].split('E')[0]),1,r_kpc)
r_kpc_2 = np.apply_along_axis(lambda a: (a[0].split('E')[1]),1,r_kpc)

#function to calculate the value
def Calculator(r_1, r_2, delta_r_1, delta_r_2, V_1, V_2, delta_V_1, delta_V_2):
    
    #Making empty arrays
    DELTA_V_KMS = np.array([])
    V_KMS = np.array([])
    R_KPC = np.array([])
    DELTA_R_KPC = np.array([])
    
    #Converting the inputs to floats
    r_1 = r_1.astype(float)
    r_2 = r_2.astype(float)
    delta_r_1 = delta_r_1.astype(float)
    delta_r_2 = delta_r_2.astype(float)
    V_1 = V_1.astype(float)
    V_2 = V_2.astype(float)
    delta_V_1 = delta_V_1.astype(float)
    delta_V_2 = delta_V_2.astype(float)
    
    #Calculating the actual value of the given data
    for a in range(0, w-(n+1)):
        Value = r_1[a]*10**r_2[a]
        R_KPC = np.append(R_KPC, Value)

    for a in range(0, w-(n+1)):
        Value = delta_r_1[a]*10**delta_r_2[a]
        DELTA_R_KPC = np.append(DELTA_R_KPC, Value)    
 
    for a in range(0, w-(n+1)):
        Value = V_1[a]*10**V_2[a]
        V_KMS = np.append(V_KMS, Value)
    
    for a in range(0, w-(n+1)):
        Value = delta_V_1[a]*10**delta_V_2[a]
        DELTA_V_KMS = np.append(DELTA_V_KMS, Value)
    return R_KPC, DELTA_R_KPC, V_KMS, DELTA_V_KMS

R_KPC, DELTA_R_KPC, V_KMS, DELTA_V_KMS = Calculator(r_kpc_1, r_kpc_2, delta_r_kpc_1, delta_r_kpc_2, V_kms_1, V_kms_2, delta_V_kms_1, delta_V_kms_2)

#TEST VALUES
G_kpc3_gy2sm = 4.49*10**(-6)

#THIS FUNCTION TAKES IN THE EXPECTED a,M,and R FOR BULGE,DISC,HALO AND RETURNS VCI
def potential_calc(R,M_i,a_i):
    G = G_kpc3_gy2sm
    return((R**2)*G*M_i)/np.sqrt((R**2+a_i**2)**3)

#EXPECTED VALUES, BULGE 0-3KPCS
M_b_sm = 1.5*10**10
a_b_kpc = 0.4

#EXPECTED VALUES, DISC 3-30KPCS
M_d_sm = 6*10**10
a_d_kpc = 5

#EXPECTED VALUES, HALO 30-40 KPCS
M_h_sm = 1.3*10**11
a_h_kpc = 12

#CHI SQUARED FUNCTION TAKES IN EXPECTED VALUES, RUNS THEM THROUGH POTENTIAL CALC, ADDS THEM IN QUADRATURE
def chisqfunc(params):
    chisq = 0
    params = [M_b_sm,a_b_kpc,M_d_sm,a_d_kpc,M_h_sm,a_h_kpc]
    #ITERATE THROUGH ENTIRE DATA SET
    for i in range(0,len(V_KMS)):
        #CALCULATE VC_tot USING potential_calc FUNCTION ABOVE
        Vc_tot=np.sqrt(potential_calc(R_KPC[i],M_b_sm,a_b_kpc)+potential_calc(R_KPC[i],M_d_sm,a_d_kpc)+potential_calc(R_KPC[i],M_h_sm,a_h_kpc))
        #CALCULATE CHISQ VALUE
        chisq += (V_KMS[i]-Vc_tot)**2/(DELTA_V_KMS[i]**2)
    return chisq

#APPENDING VALUES TO ARRAY FOR OPTIMIZATION
x0 = np.array([M_b_sm, M_d_sm, M_h_sm, a_b_kpc, a_d_kpc, a_h_kpc])

#OPTIMIZATION USING SCIPY.minimize
optimize = optimize.minimize(chisqfunc, x0,method='Nelder-Mead')
print(optimize.x)
#EXTRACTING OPTIMIZED VALUES
M_b_sm1 = optimize.x[0]
a_b_kpc1 = optimize.x[3]
M_d_sm1 = optimize.x[1]
a_d_kpc1 = optimize.x[4]
M_h_sm1 = optimize.x[2]
a_h_kpc1 = optimize.x[5]

#ADDING VALUES IN QUADRATURE, TAKES IN OPTIMIZED PARAMETERS
def add_quad(optimize):
    V_tot = []
    #ITERATING THROUGH ENTIRE RANGE OF VALUES
    for i in range(0,len(V_KMS)):
        #ADDING IN QUADRATURE
        V_C = np.sqrt(potential_calc(R_KPC[i],M_b_sm1,a_b_kpc1) + potential_calc(R_KPC[i],M_d_sm1,a_d_kpc1) + potential_calc(R_KPC[i],M_h_sm1,a_h_kpc1))
        #APPENDING TO LIST
        V_tot.append(V_C)
    return V_tot  

#PLOTTING
fig1 = plt.plot(R_KPC,V_KMS)
plt.plot(R_KPC,add_quad(optimize))
plt.title("$V_c$ vs R [Milky Way Galaxy]",fontsize = 16)
plt.xlabel(r"$Radius \quad [kpc]$")
plt.ylabel(r"$V_c \quad [\frac{km}{s}]$")
plt.legend(loc='lower center',fontsize=10)      
plt.show()

#PRINT FUNCTION
print('Figure 1: The rotation curve for the Milky Way Galaxy using observed values and our model. My model determined the following values: M_bulge = ',M_b_sm1, 'Solar Masses, M_disk =', M_d_sm1, 'Solar Masses, M_halo =', M_h_sm1, 'Solar Masses, a_bulge =', a_b_kpc1, 'kpc, a_disk =', a_d_kpc1, 'kpc, and a_halo =', a_h_kpc1, 'kpc')
