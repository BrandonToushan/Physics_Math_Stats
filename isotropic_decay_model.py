####MODELING ISOTROPIC DECAY OF Mn-54

#The dominant form of decay of Mn-54 is electron capture, followed by gamma partical emission
#to the 834.855 keV excited state of Cr-54 with a branching ratio of  99.9977% of the time this occuring
#while 0.0003% of the time weak electron transfer occurs with beta plus transfer

#The attenuation length for this decay is 27.5 microns

# %. http://www.nucleide.org/DDEP_WG/Nuclides/Mn-54_tables.pdf
#https://physics.nist.gov/cgi-bin/ffast/ffast.pl?Formula=MnO2&gtype=5&lower=&upper=&density=0.01 Attenuation length.. 27.5 microns

#http://web-docs.gsi.de/~stoe_exp/web_programs/x_ray_absorption/index.php

import math 

from numpy import random

import numpy as np

from matplotlib import pyplot

from mpl_toolkits.mplot3d import Axes3D

##DEFINE A BASIC CLASS FOR SOURCE VECTOR MODELLING

class aSource(object):
    
    #CONSTRUCTOR OF EMPTY SOURCE CLASS
    
    def __init__(self, n=1, title = "Empty Source", orig=[0,0,0]):
        self.nDraws = n
        self.title = title
        self.origin = orig
        self.draw()
        
    #GENERATE VECTORS (IN THIS CASE NONE)
        
    def draw(self):
        self.x = None
        self.y = None
        self.v = None
        self.rand_angles = None

    #THIS METHOD RETURNS THE DATA FROM THE CLASS
        
    def get(self):
        return {"Title": self.title, "Vectors": self.v, "Xs": self.x, "Ys": self.y, "Angles (rad)": self.rand_angles, "Origin": self.origin}
    
#TEST EMPTY CLASS

#print "--------------Empty Class Test---------------------"

a_Source = aSource()

#print a_Source.get()

##CREATE ISOTROPIC DECAY VECTORS, USING INHERITATNCE FROM ABOVE aSOURCE CLASS

class isoSource(aSource):

    def __init__(self, n = 1, orig = [0,0,0]):
        aSource.__init__(self, title = "2D Isotropic Source", n = n, orig = orig)

    #OVERRIDES EMPTY VECTOR GENERATION FROM aSOURCE
        
    def draw(self):
        self.rand_angles = random.uniform(0,2*np.pi,self.nDraws)
        self.x = [math.cos(i) for i in self.rand_angles]
        self.y = [math.sin(i) for i in self.rand_angles]
        self.v = [[x,y] for x,y in zip(self.x, self.y)]

#PRINT ISOTROPIC SOURCE VECTORS
        
#print "--------------Isotropic 2D-Source-------------------"
        
iso_source = isoSource()

#print iso_source.get()

##DEFINE A CLASS FOR PLANE SHAPED ORIGIN OF ISOTROPIC DECAY VECTORS ABOVE

class planeSource(object):

        #NUMBER OF SOURCES PRESENT & SOURCE DIMENSIONS
    
	def __init__(self, numSource = 1, x_max = 1, y_max = 1, z_max = 1):
            self.numSources = numSource
            self.x_max = x_max
            self.y_max = y_max
            self.z_max = z_max
            self.draw()
	
	#USING ISOTROPIC DECAY SOURCE CLASS FROM ABOVE
		
	def draw(self):
            self.isoSor=[isoSource(n=1, orig=[random.uniform(0,self.x_max),random.uniform(0,self.y_max),random.uniform(0,self.y_max)]) for i in range(0, self.numSources)]

            #DETERMINE ORIGIN AND VECTORS FOR ISOTROPIC DECAY
            self.origins=[sor.origin for sor in self.isoSor]
            self.vectors=[sor.v for sor in self.isoSor]

#PRINT OUTPUTS
		
#print "--------------Isotropic Decay Vectors from Plane Shaped Source-------------------"

final_plane = planeSource()

#Radiation Source limits
xmin = 0
xmax = 1
ymin = 0
ymax = 1
zmin = 0
zmax = 0.05

def source_exit (xmin,xmax,ymin,ymax,zmin,zmax,vec_origin,vec_direct):
     n_zbot = (zmin-vec_origin[2])/vec_direct[2]
     n_ztop = (zmax-vec_origin[2])/vec_direct[2]
     n_xmin = (xmin-vec_origin[0])/vec_direct[0]
     n_xmax = (xmax-vec_origin[0])/vec_direct[0]
     n_ymin = (ymin-vec_origin[1])/vec_direct[1]
     n_ymax = (ymax-vec_origin[1])/vec_direct[1]
     mylist = [n_zbot,n_ztop,n_xmin,n_xmax,n_ymin,n_ymax]
     
     low_n = min(i for i in mylist if i >= 0)
     #exit_plane = mylist.index(low_n)
     full_vector = [i * low_n for i in vec_direct]
     #vectorpnt = [vec_origin[0] + full_vector[0],vec_origin[1] + full_vector[1], vec_origin[2] + full_vector[2]]
     d_travelled = np.linalg.norm(full_vector)
     return d_travelled

distances = []

for x in range(1000):
    final_plane = planeSource()
    zvec = np.cos(random.uniform(0,2*np.pi))
    Vector_Direction = [final_plane.vectors [0][0][0],final_plane.vectors [0][0][1],zvec]
    Origin_vector = [final_plane.origins [0][0] * xmax,final_plane.origins [0][1] * ymax,final_plane.origins[0][2] * zmax]
    y = source_exit(xmin,xmax,ymin,ymax,zmin,zmax,Origin_vector,Vector_Direction)
    distances.append(y)
figa = pyplot.hist(distances, bins=100)
pyplot.title("Distribution of Lengths in Source" )
pyplot.xlabel('Position in Source')
pyplot.ylabel('Frequency')
pyplot.show()

#REUSE PLANESOURCE FROM ABOVE TO RUN FUNCTION AGAIN WITH 10000 POINTS FOR PLOTTING, NECCESSARY DUE TO ISSUES COMBINING CODE BETWEEN PARTNERS

class planeSource_1(object):

        #NUMBER OF SOURCES PRESENT & SOURCE DIMENSIONS
    
	def __init__(self, numSource_a = 10000, x_max = 1, y_max = 1, z_max = 1):
            self.numSources_a = numSource_a
            self.x_max = x_max
            self.y_max = y_max
            self.z_max = z_max
            self.draw()
	
	#USING ISOTROPIC DECAY SOURCE CLASS FROM ABOVE
		
	def draw(self):
            self.isoSor_a=[isoSource(n = 1, orig=[random.uniform(0,self.x_max),random.uniform(0,self.y_max),random.uniform(0,self.y_max)]) for i in range(0, self.numSources_a)]

            #DETERMINE ORIGIN AND VECTORS FOR ISOTROPIC DECAY
            
            self.origins = [sor.origin for sor in self.isoSor_a]
            self.vectors = [sor.v for sor in self.isoSor_a]


final_plane_a = planeSource_1()
            
vec_origin_a = final_plane_a.origins
vec_direct_a = final_plane_a.vectors

#print vec_origin, vec direct

#EXTRACT X,Y,Z COARDINATES FROM SOURCE DISTRIBUTION FOR PLOTTING

x_a = [item[0] for item in vec_origin_a]
y_a = [item[1] for item in vec_origin_a]
z_a = [item[2] for item in vec_origin_a]


fig1 = pyplot.figure()
ax = Axes3D(fig1)
ax.set_title("Positions of Isotropic Decay Sources (x,y,z)")
ax.set_xlabel("x Coordinate (cm)")
ax.set_ylabel("y Coordinate (cm)")
ax.set_zlabel("z Coordinate (cm)")
ax.scatter(x_a,y_a,z_a)
pyplot.show()

##PLOTS OF ORIGIN POSITIONS OF 3D SOURCES

#x-y COORDINATE PLOT

fig2 = pyplot.scatter(x_a,y_a,s=5)
pyplot.title("Positions of Isotropic Decay Sources (x,y)")
pyplot.xlabel('x Coordinates (cm)')
pyplot.ylabel('y Coordinates (cm)')
pyplot.show()

#y-z COORDINATE PLOT

fig3 = pyplot.scatter(y_a,z_a,s=5)
pyplot.title("Positions of Isotropic Decay Sources (y,z)")
pyplot.xlabel('y Coordinates (cm)')
pyplot.ylabel('z Coordinates (cm)')
pyplot.show()

#z-x COORDINATE PLOT

fig4 = pyplot.scatter(z_a,x_a,s=5)
pyplot.title('Positions of Isotropic Decay Sources (z,x)')
pyplot.xlabel('z Coordinates (cm)')
pyplot.ylabel('x Coordinates (cm)')
pyplot.show()

#FRACTION OF X-RAYS THAT ESCAPE

#VALUES
att_l_m = (27.5)*10**(-6)
h_m = (0.01)*10**(-6)

#CREATE FUNCTION TO FIND AMOUNT OF PHOTONS ESCAPING
def escape(self,att_l_m):

    #COUNTER FOR ESCAPE
    a = 0

    #RANDOM DISTRIBUTION OF ESCAPES
    dist = [random.uniform(0,1) for i in range(0,self.n)]

    #IF PHOTON EXCEEDS BOUNDARY, ADD 1 TO COUNT, RETURN RATIO
    for i in range(self.n):

            if dist[i] < np.exp(-inter[i]/att_l_m):
                       a += 1
