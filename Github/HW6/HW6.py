#!/usr/bin/env python
# coding: utf-8

# In[46]:


# from numpy import math
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import matplotlib
import math
import pandas as pd
np.core
from sys import maxsize
from numpy import set_printoptions
set_printoptions(threshold=maxsize)

print("hi")





#*****************Constants********************#
Lx = 2000
Nx = 101
Ly = 2000
Ny = 101
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
print("dx = ",dx)
print("dy = ",dy)


#r = 500.0
tol = 10**-7

a = 180
b = 600

x_int_1 = (Lx/2) + (b/2)
y_int_1 = Ly/2
x_int_2 = (Lx/2) - (b/2)
y_int_2 = Ly/2

Ttot = 20 #DAYS

#****************************************#
#*************Functions******************#
#****************************************#

def forceEqual(array): #left right symmetry
    finarr = np.transpose(array)
    i = 0
    while (i != 52):
        array[i] = array[102-i]
        i+=1
    return np.transpose(finarr)

def twoDSum (array):
    newarr = []
    for i in array:
        newarr+=[math.fsum(i)]
    return (math.fsum(newarr))


def U_func (PHI):
    U = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            U[i][j] = (PHI[i][j-1]-PHI[i][j+1])/(2*dy)
    return (U)
    
def V_func (PHI):
    V = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            V[i][j] = (PHI[i+1][j]-PHI[i-1][j])/(2*dx)
    return (V)


def J_1_func (PHI, VORT): #return standard size Nx+2 x Ny+2 array
    J_1 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1, Nx+1):
        for j in range(1,Ny+1):
            J_1[i][j] = (1/(2*dx))*(PHI[i+1][j]-PHI[i-1][j])*(1/(2*dy))*(VORT[i][j+1] - VORT[i][j-1]) - (1/(2*dy))*(PHI[i][j+1] - PHI[i][j-1])*(1/(2*dx))*(VORT[i+1][j]-VORT[i-1][j])
    return (J_1)

def J_2_func (PHI, VORT):
    J_2 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_2[i][j] = (1/(2*dx))*(PHI[i+1][j]*(1/(2*dy))*(VORT[i+1][j+1]-VORT[i+1][j-1]) - PHI[i-1][j]*(1/(2*dy))*(VORT[i-1][j+1]-VORT[i-1][j-1])) - (1/(2*dy))*(PHI[i][j+1]*(1/(2*dx))*(VORT[i+1][j+1]-VORT[i-1][j+1]) - PHI[i][j-1]*(1/(2*dx))*(VORT[i+1][j-1]-VORT[i-1][j-1]))    
    return (J_2)


def J_3_func (PHI, VORT):
    J_3 =  np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_3[i][j] = (1/(2*dy))*(VORT[i][j+1]*(1/(2*dx))*(PHI[i+1][j+1]-PHI[i-1][j+1]) - VORT[i][j-1]*(1/(2*dx))*(PHI[i+1][j-1] - PHI[i-1][j-1])) - (1/(2*dx))*(VORT[i+1][j]*(1/(2*dy))*(PHI[i+1][j+1]-PHI[i+1][j-1]) - VORT[i-1][j]*(1/(2*dy))*(PHI[i-1][j+1]-PHI[i-1][j-1]))
    return (J_3)

def MaxCourant(PHI, dt, dx ): #dx = 20
    absmaxU = np.amax(np.absolute(U_func(PHI)[2:-2,2:-2]))
    absmaxV = np.amax(np.absolute(V_func(PHI)[2:-2,2:-2]))
    return (absmaxU*dt/dx + absmaxV*dt/dx)

#IMPORTANT, RANGE (ghost nodes + walls ignored)
def totalKE (PHI):
    total = 0
    U = np.copy(U_func(PHI))
    V = np.copy(V_func(PHI))
    for i in range(2,Nx):
        for j in range(2,Ny):
            total += 0.5*(U[i][j]**2 + V[i][j]**2)
    return (total)
#IMPORTANT, VORTICITY AT WALL???????
def enstrophy (VORT):
    total = 0
    for i in range(2,101):
        for j in range(2,101):
            total += (VORT[i][j])**2
    return (total)

#*************Vorticity Initilaize******************#

VORT = np.zeros((Nx+2, Ny+2), dtype = float)

circ_factor = 8*10**(-5)

for i in range(1,Nx+1): #101 physical node
    for j in range(1,Ny+1):
        VORT[i][j] = circ_factor*math.exp(-(((i-1)*dx - x_int_1)**2 + ((j-1)*dy-y_int_1)**2)/a**2) + circ_factor*math.exp(-(((i-1)*dx - x_int_2)**2 + ((j-1)*dy-y_int_2)**2)/a**2)

#Boundary Condition
for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = VORT[i][2]
    VORT[i][0] = VORT[i][Ny-1]
    VORT[i][1] = VORT[i][Ny]
print("finish VORT initialize")



#print(VORT[1:-1,1:-1])


plt.figure(0)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.savefig("INT,UNNORM-VORT-run_B.svg")



VORT_Avg = twoDSum(VORT[1:-1,1:-1])/(VORT[1:-1,1:-1].size)
print("mean vorticity = ", VORT_Avg)

"""
VORT_loop_counter = 0

        #AVERGING
while (VORT_Avg >  tol):
    for i in range(0,Nx+2):
        for j in range(0,Ny+2):
            VORT[i][j] = VORT[i][j] - VORT_Avg
            
    VORT_Avg = twoDSum(VORT[1:-1,1:-1])/(VORT[1:-1,1:-1].size)
    VORT_loop_counter+=1

print("VORT loop counter", VORT_loop_counter)
""" 


plt.figure(1)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
#plt.savefig("INT,NORM-VORT-run_B.svg")
#VORT = np.transpose(VORT_R)



#*****************Advancing Lists********************#

PHI = np.zeros((Nx+2, Ny+2), dtype = float)
Residual = np.zeros((Nx+2, Ny+2), dtype = float)


U = np.zeros((Nx+2, Ny+2), dtype = float)
V = np.zeros((Nx+2, Ny+2), dtype = float)
#Initial residual:             
for i in range( 0,Nx+2): 
    for j in range(0,Ny+2):
        Residual[i][j] =  -(VORT[i][j])

plt.figure(100)
plt.imshow( np.transpose(Residual[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")




"""
for a in range(0,Ny+2):
    Residual[Nx+1][a] = Residual[2][a]
    Residual[0][a] = Residual[Nx-1][a]
    Residual[1][a] = Residual[Nx][a]
for b in range(0,Nx+2):
    Residual[b][Ny+1] = Residual[b][2]
    Residual[b][0] = Residual[b][Ny-1]
    Residual[b][1] = Residual[b][Ny]

"""    
  
print("finish Residual initialization")


#test #print(u_cur.shape)
#result = 53x53

###########################Criterion for Convergence#######################

EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
print("This is initial EPS = " ,EPS)

########################### ALPHA #######################   

#beta = dx/dy
#F-PLANE APPROXIMATION
sigma = (1/(1+(dx/dy)**2))*(math.cos(math.pi/Nx) + math.cos(math.pi/Ny)*(dx/dy)**2)
alpha = 2/(1+math.sqrt(1-sigma**2))

print("sigma =", sigma)    
print("alpha =", alpha)

timeStep = []
EPScounter = []

#***************MAIN LOOP*********************#
maincounter = 0
print("entering main loop")

while (EPS > tol):

    for j in range(2,Ny):  
    
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        
        for i in range(2,Nx): #101 physical node WALL
            
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finish

            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    
    maincounter +=1
    timeStep += [maincounter]
    EPScounter +=[EPS]

#URN
PHI = forceEqual(PHI)
PHI = np.transpose(forceEqual(np.transpose(PHI)))

        
#print(EPS)        
  
 
    
#print((PHI[1:-1,1:-1]))
    
    


#NUMERICAL VORTICITY
VORT_NUM = np.zeros((Nx+2, Ny+2), dtype = float)
for i in range(1, Nx+1):
    for j in range(1,Ny+1):
        VORT_NUM[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j])+(1/dy**2)*(PHI[i][j-1]-2*PHI[i][j]+ PHI[i][j+1])
        
plt.figure(5)


x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
                            #URN
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,2000,0,2000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)


#quiver plot


U = np.copy(U_func(PHI))
V = np.copy(V_func(PHI))


plt.figure(4)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(PHI[1:-1,1:-1]), interpolation='none', extent=[0,2000,0,2000])
plt.colorbar()#URN
plt.quiver(X, Y, np.transpose(U[1:-1,1:-1]), np.transpose(V[1:-1,1:-1]))        
plt.savefig('quiverB.svg', dpi = 300)



plt.figure(10)


VORT_DIFF = np.zeros((Nx+2, Ny+2), dtype = float)
for i in range(1,Nx+1):
    for j in range(1,Ny+1):
        VORT_DIFF[i][j]=VORT_NUM[i][j]-VORT[i][j]
        
plt.imshow( np.transpose(np.absolute(VORT_DIFF[1:-1,1:-1])), interpolation='none', extent=[0,2000,0,2000])
#plt.colorbar(ticks = np.linspace(0,4.5,18,endpoint=True))
plt.colorbar()
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
#plt.savefig('numdiff_RUNB.svg', dpi = 300)

#print(VORT_DIFF)

plt.figure(11)
plt.plot(timeStep, EPScounter)
plt.xscale("linear")
plt.xlabel("Time Step")
plt.ylabel("EPS")
plt.yscale("log")
plt.title("Epsilon - Time Step Graph [RUN B]")

#plt.savefig("EPS-TIME_RUNB.svg")

print("Total KE =", totalKE(PHI))



# In[94]:


# from numpy import math
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import matplotlib
import math
import pandas as pd
np.core
from sys import maxsize
from numpy import set_printoptions
set_printoptions(threshold=maxsize)

print("hi")





#*****************Constants********************#
Lx = 2000*1000
Nx = 101
Ly = 2000*1000
Ny = 101
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
print("dx = ",dx)
print("dy = ",dy)


#r = 500.0
tol = 10**-7

a = 180*1000
b = 600*1000

x_int_1 = (Lx/2) + (b/2)
y_int_1 = Ly/2
x_int_2 = (Lx/2) - (b/2)
y_int_2 = Ly/2

Ttot = 1728000 #SECONDS
dt = 1000
Nt = Ttot/dt




#****************************************#
#*************Functions******************#
#****************************************#

def forceEqual(array): #left right symmetry
    finarr = np.transpose(array)
    i = 0
    while (i != 52):
        array[i] = array[102-i]
        i+=1
    return np.transpose(finarr)

def twoDSum (array):
    newarr = []
    for i in array:
        newarr+=[math.fsum(i)]
    return (math.fsum(newarr))


def U_func (PHI):
    dy = 20000
    U = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            U[i][j] = (PHI[i][j-1]-PHI[i][j+1])/(2*dy)
    return (U)
    
def V_func (PHI):
    dx = 20000
    V = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            V[i][j] = (PHI[i+1][j]-PHI[i-1][j])/(2*dx)
    return (V)


def J_1_func (PHI, VORT): #return standard size Nx+2 x Ny+2 array
    J_1 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1, Nx+1):
        for j in range(1,Ny+1):
            J_1[i][j] = (1/(2*dx))*(PHI[i+1][j]-PHI[i-1][j])*(1/(2*dy))*(VORT[i][j+1] - VORT[i][j-1]) - (1/(2*dy))*(PHI[i][j+1] - PHI[i][j-1])*(1/(2*dx))*(VORT[i+1][j]-VORT[i-1][j])
    return (J_1)

def J_2_func (PHI, VORT):
    J_2 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_2[i][j] = (1/(2*dx))*(PHI[i+1][j]*(1/(2*dy))*(VORT[i+1][j+1]-VORT[i+1][j-1]) - PHI[i-1][j]*(1/(2*dy))*(VORT[i-1][j+1]-VORT[i-1][j-1])) - (1/(2*dy))*(PHI[i][j+1]*(1/(2*dx))*(VORT[i+1][j+1]-VORT[i-1][j+1]) - PHI[i][j-1]*(1/(2*dx))*(VORT[i+1][j-1]-VORT[i-1][j-1]))    
    return (J_2)


def J_3_func (PHI, VORT):
    J_3 =  np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_3[i][j] = (1/(2*dy))*(VORT[i][j+1]*(1/(2*dx))*(PHI[i+1][j+1]-PHI[i-1][j+1]) - VORT[i][j-1]*(1/(2*dx))*(PHI[i+1][j-1] - PHI[i-1][j-1])) - (1/(2*dx))*(VORT[i+1][j]*(1/(2*dy))*(PHI[i+1][j+1]-PHI[i+1][j-1]) - VORT[i-1][j]*(1/(2*dy))*(PHI[i-1][j+1]-PHI[i-1][j-1]))
    return (J_3)

def J_final (PHI,VORT):
    return ((1/3)*(J_3_func(PHI,VORT)+J_2_func(PHI,VORT)+J_1_func(PHI,VORT)))



def MaxCourant(PHI, dt,dx ): #dx = 20
    absmaxU = np.amax(np.absolute(U_func(PHI)[2:-2,2:-2]))
    print("Abs max U is", absmaxU)
    absmaxV = np.amax(np.absolute(V_func(PHI)[2:-2,2:-2]))
    print("Abs max V is", absmaxV)
    return (absmaxU*dt/dx + absmaxV*dt/dx)

#IMPORTANT, RANGE (ghost nodes + walls ignored)
def totalKE (PHI):
    total = 0
    U = np.copy(U_func(PHI))
    V = np.copy(V_func(PHI))
    for i in range(2,Nx):
        for j in range(2,Ny):
            total += 0.5*(U[i][j]**2 + V[i][j]**2)
    return (total)
#IMPORTANT, VORTICITY AT WALL???????
def enstrophy (VORT):
    total = 0
    for i in range(2,101):
        for j in range(2,101):
            total += (VORT[i][j])**2
    return (total)
#NUMERICAL VORTICITY
def VORT_NUM_func (PHI):
    VORT_NUM = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1, Nx+1):
        for j in range(1,Ny+1):
            VORT_NUM[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j])+(1/dy**2)*(PHI[i][j-1]-2*PHI[i][j]+ PHI[i][j+1])
    return (VORT_NUM)
        

#******************************************************************************#
#***********************Initialize Vorticity***********************************#
#******************************************************************************#


VORT = np.zeros((Nx+2, Ny+2), dtype = float)

circ_factor = 8*10**(-5)

for i in range(1,Nx+1): #101 physical node
    for j in range(1,Ny+1):
        VORT[i][j] = circ_factor*math.exp(-(((i-1)*dx - x_int_1)**2 + ((j-1)*dy-y_int_1)**2)/a**2) + circ_factor*math.exp(-(((i-1)*dx - x_int_2)**2 + ((j-1)*dy-y_int_2)**2)/a**2)

        
#Boundary Condition

for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    VORT[i][1] = VORT[i][Ny]




print("finish VORT initialize")

#print(VORT[1:-1,1:-1])


plt.figure(0)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.title("Initial Vorticity")
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
#plt.savefig("INT,UNNORM-VORT-run_B.svg")






#******************************************************************************#
#*********************************** NORMALIZATION ****************************#
#******************************************************************************#

VORT_Avg = twoDSum(VORT[1:-1,1:-1])/(VORT[1:-1,1:-1].size)
VORT_loop_counter = 0

print("VORT loop counter", VORT_loop_counter)
print("mean vorticity = ", VORT_Avg)

"""
plt.figure(1)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
"""
#plt.savefig("INT,NORM-VORT-run_B.svg")
#VORT = np.transpose(VORT_R)


#******************************************************************************#
#*****************POISSON ADVANCING LISTS**************************************#
#******************************************************************************#

PHI = np.zeros((Nx+2, Ny+2), dtype = float)
Residual = np.zeros((Nx+2, Ny+2), dtype = float)


U = np.zeros((Nx+2, Ny+2), dtype = float)
V = np.zeros((Nx+2, Ny+2), dtype = float)
#Initial residual:             
for i in range( 0,Nx+2): 
    for j in range(0,Ny+2):
        Residual[i][j] =  -(VORT[i][j])

plt.figure(100)
plt.imshow( np.transpose(Residual[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()
plt.title("Initial Residual")
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")




"""
for a in range(0,Ny+2):
    Residual[Nx+1][a] = Residual[2][a]
    Residual[0][a] = Residual[Nx-1][a]
    Residual[1][a] = Residual[Nx][a]
for b in range(0,Nx+2):
    Residual[b][Ny+1] = Residual[b][2]
    Residual[b][0] = Residual[b][Ny-1]
    Residual[b][1] = Residual[b][Ny]

"""    
  
print("finish Residual initialization")


#test #print(u_cur.shape)
#result = 53x53
#*************************************************************************#
###########################Criterion for Convergence#######################
#*************************************************************************#
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
print("This is initial EPS = " ,EPS)


#******************************************************************************#
#************************************ALPHA*************************************#
#******************************************************************************#

#beta = dx/dy
#F-PLANE APPROXIMATION
sigma = (1/(1+(dx/dy)**2))*(math.cos(math.pi/Nx) + math.cos(math.pi/Ny)*(dx/dy)**2)
alpha = 2/(1+math.sqrt(1-sigma**2))

print("sigma =", sigma)    
print("alpha =", alpha)

timeStep = []
EPScounter = []

#******************************************************************************#
#***********************ADVANCING LIST (JACOBIANS)*****************************#
#******************************************************************************#
J_older = np.zeros((Nx+2, Ny+2), dtype = float) #n-2
J_old = np.zeros((Nx+2, Ny+2), dtype = float)   #n-1
J_cur = np.zeros((Nx+2, Ny+2), dtype = float)   #n
J_new = np.zeros((Nx+2, Ny+2), dtype = float)   #n+1

VORT_new = np.zeros((Nx+2, Ny+2), dtype = float)



#******************************************************************************#
#*********************** 0th Solver        ************************************#
#******************************************************************************#

EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0
maincounter = 0
while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    
    maincounter +=1
    #timeStep += [maincounter]
    #EPScounter +=[EPS]

print("0th Step: SOR Iteration Counter =", maincounter )



U = np.copy(U_func(PHI))
V = np.copy(V_func(PHI))


plt.figure(10000)
plt.title("Velocity Plot original")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(U[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()  

VORT_NUM = np.copy(VORT_NUM_func(PHI))
plt.figure(10020)
plt.title("numerical vort init")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none',extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()  


#plt.savefig("NumVort_RUN_B.svg", dpi = 300)

#******************************************************************************#
#*****************K total, Ensthropy total, max Courant************************#
#******************************************************************************#


print("0th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))
#print("0th Enstropy =", enstrophy(VORT))
      
#******************************************************************************#
#*****************************FIRST TIME STEP**********************************#
#******************************************************************************#

#Jacobian
J_older = np.copy(J_final(PHI,VORT))

#Step in time
for i in range(0,Nx+2):
    for j in range(0,Ny+2):
        VORT[i][j] = VORT[i][j] - J_older[i][j]*dt

#Boundary Condition
for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    #VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    #VORT[i][1] = VORT[i][Ny]

plt.figure(50)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
plt.title("I'm VORT")                     #URN
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()

    
#####POISSON SOLVER######

#ORDER is ******* IMPORTANT
PHI = np.zeros((Nx+2, Ny+2), dtype = float)
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0


maincounter = 0

while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    maincounter +=1
    
print("1th Step: SOR Iteration Counter =", maincounter )        
 
#UPDATING QUANTITIES
print("First Iteration Results:::::")
print("1th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))
#print("1st Enstropy =", enstrophy(VORT))
     



plt.figure(1493)

plt.title("I'm VORT_NUM after 1st EF")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)

VORT_NUM = np.copy(VORT_NUM_func(PHI))  
#URN
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)




#******************************************************************************#
#*****************************SECOND TIME STEP*********************************#
#******************************************************************************#


#Jacobian
J_old = np.copy(J_final(PHI,VORT))

#Step in time
for i in range(0,Nx+2):
    for j in range(0,Ny+2):
        VORT[i][j] = VORT[i][j] - J_old[i][j]*dt

#Boundary Condition
for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    #VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    #VORT[i][1] = VORT[i][Ny]

plt.figure(5012)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
plt.title("I'm VORT (2nd EF)")                     #URN
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()

    
#####POISSON SOLVER######

#ORDER is ******* IMPORTANT
PHI = np.zeros((Nx+2, Ny+2), dtype = float)
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0

maincounter = 0

while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    maincounter +=1
    
print("1th Step: SOR Iteration Counter =", maincounter )        
 
#UPDATING QUANTITIES
print("First Iteration Results:::::")
print("2th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))

#print("2nd Enstropy =", enstrophy(VORT))


#print("Enstropy =", enstrophy(VORT_NUM(PHI)))       



plt.figure(14931)

plt.title("I'm VORT_NUM after 2st EF")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)

VORT_NUM = np.copy(VORT_NUM_func(PHI))  
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)


#******************************************************************************#
#****************************** MAIN TIME LOOP ********************************#
#******************************************************************************#

remove_counter=0

while(remove_counter < 300):

    #Jacobian
    J_old = np.copy(J_final(PHI,VORT))

    #Step in time
    for i in range(0,Nx+2):
        for j in range(0,Ny+2):
            VORT[i][j] = VORT[i][j] - J_old[i][j]*dt

    #Boundary Condition
    for j in range(0,Ny+2):
        VORT[Nx+1][j] = 0
        VORT[0][j] = 0
        #VORT[1][j] = VORT[Nx][j]
    for i in range(0,Nx+2):
        VORT[i][Ny+1] = 0
        VORT[i][0] = 0






    #####POISSON SOLVER######

    #ORDER is ******* IMPORTANT
    PHI = np.zeros((Nx+2, Ny+2), dtype = float)
    EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    # = 1.0

    maincounter = 0

    while (EPS > tol):
        for j in range(2,Ny):  
            #update ghost nodes along phi[i][Ny+1]
            if (j == Ny):  #IMPORTANT , NY-1 for wall
                for z in range(0,Nx+2): 
                    PHI[z][Ny+1] = PHI[z][2]
            #update finish
            for i in range(2,Nx): #101 physical node WALL
                #update ghost nodes along phi[Nx+1][j]
                if (i == Nx):
                    for k in range(0,Ny+2):
                        PHI[Nx+1][k] = PHI[2][k]
                #update finisH
                Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
                PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))

            #Update ghosts at PHI[0][j]
            for k in range(0,Ny+2):
                PHI[0][k] = PHI[Nx-1][k]
                PHI[1][k] = PHI[Nx][k] #make sure left right equal

        #update finish
        for t in range(0,Nx+2):
            PHI[t][0] = PHI[t][Ny-1]
            PHI[t][1] = PHI[t][Ny] 
        #make sure left right equal

        #update finish
        #Update ghosts at PHI[i][0]
        tempResidual = np.copy(Residual)
        tempResidual = np.absolute(tempResidual[2:-2,2:-2])

        EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
        #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
        #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))

        maincounter +=1

    print("big counter =", remove_counter)
    remove_counter+=1
    


    if (remove_counter%10 == 0):
        
        plt.figure(remove_counter+102432)
        x = np.linspace(0,2000,101)
        y = np.linspace(0,2000,101)
        plt.xlabel("X Domain [m]")
        plt.ylabel("Y Domain [m]")
        X,Y = np.meshgrid(x,y)
        plt.title("I'm VORT (final nd EF)")                     #URN
        plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
        plt.colorbar()



        plt.figure(remove_counter+102433)
        x = np.linspace(0,2000,101)
        y = np.linspace(0,2000,101)
        plt.xlabel("X Domain [m]")
        plt.ylabel("Y Domain [m]")
        X,Y = np.meshgrid(x,y)
        plt.title("I'm PSI (final nd EF)")                     #URN
        plt.imshow( np.transpose(PHI[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
        plt.colorbar()
        print("remove_counter =", remove_counter)
        print("nth Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
        print("Total KE =", totalKE(PHI))













"""
maincounter = 0
print("entering main loop")

while (EPS > tol):

    for j in range(2,Ny):  
    
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        
        for i in range(2,Nx): #101 physical node WALL
            
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finish

            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    
    maincounter +=1
    timeStep += [maincounter]
    EPScounter +=[EPS]

#URN
#PHI = forceEqual(PHI)
#PHI = np.transpose(forceEqual(np.transpose(PHI)))

        
#print(EPS)        
  
 
    
#print((PHI[1:-1,1:-1]))
    
    
#******************************************************************************#
#*****************************NUMERICAL VORTICITY******************************#
#******************************************************************************#

VORT_NUM = np.copy(VORT_NUM(PHI))

plt.figure(5)

plt.title("numerical vorticity")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
                            #URN
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,2000,0,2000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)




#******************************************************************************#
#*********************************Velocity Plots*******************************#
#******************************************************************************#
U = np.copy(U_func(PHI))
V = np.copy(V_func(PHI))


plt.figure(4)
plt.title("Velocity Plot")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(PHI[1:-1,1:-1]), interpolation='none', extent=[0,2000,0,2000])
plt.colorbar()#URN
plt.quiver(X, Y, np.transpose(U[1:-1,1:-1]), np.transpose(V[1:-1,1:-1]))        
plt.savefig('quiverB.svg', dpi = 300)


#******************************************************************************#
#*******************************VORTICITY DIFFERENCE***************************#
#******************************************************************************#

plt.figure(10)
VORT_DIFF = np.zeros((Nx+2, Ny+2), dtype = float)
for i in range(1,Nx+1):
    for j in range(1,Ny+1):
        VORT_DIFF[i][j]=VORT_NUM[i][j]-VORT[i][j]
        
plt.imshow( np.transpose(np.absolute(VORT_DIFF[1:-1,1:-1])), interpolation='none', extent=[0,2000,0,2000])
#plt.colorbar(ticks = np.linspace(0,4.5,18,endpoint=True))
plt.colorbar()
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.title("Vorticity Difference")
#plt.savefig('numdiff_RUNB.svg', dpi = 300)

#print(VORT_DIFF)

print("Total KE =", totalKE(PHI))


print(VORT_DIFF)
"""


# In[3]:


# from numpy import math
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import matplotlib
import math
import pandas as pd
np.core
from sys import maxsize
from numpy import set_printoptions
set_printoptions(threshold=maxsize)

print("hi")





#*****************Constants********************#
Lx = 2000*1000
Nx = 101
Ly = 2000*1000
Ny = 101
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
print("dx = ",dx)
print("dy = ",dy)


#r = 500.0
tol = 10**-7

a = 180*1000
b = 600*1000

x_int_1 = (Lx/2) + (b/2)
y_int_1 = Ly/2
x_int_2 = (Lx/2) - (b/2)
y_int_2 = Ly/2

Ttot = 1728000 #SECONDS
dt = 1000
Nt = int(Ttot/dt)




#****************************************#
#*************Functions******************#
#****************************************#

def forceEqual(array): #left right symmetry
    finarr = np.transpose(array)
    i = 0
    while (i != 52):
        array[i] = array[102-i]
        i+=1
    return np.transpose(finarr)

def twoDSum (array):
    newarr = []
    for i in array:
        newarr+=[math.fsum(i)]
    return (math.fsum(newarr))


def U_func (PHI):
    dy = 20000
    U = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            U[i][j] = (PHI[i][j-1]-PHI[i][j+1])/(2*dy)
    return (U)
    
def V_func (PHI):
    dx = 20000
    V = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            V[i][j] = (PHI[i+1][j]-PHI[i-1][j])/(2*dx)
    return (V)


def J_1_func (PHI, VORT): #return standard size Nx+2 x Ny+2 array
    J_1 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1, Nx+1):
        for j in range(1,Ny+1):
            J_1[i][j] = (1/(2*dx))*(PHI[i+1][j]-PHI[i-1][j])*(1/(2*dy))*(VORT[i][j+1] - VORT[i][j-1]) - (1/(2*dy))*(PHI[i][j+1] - PHI[i][j-1])*(1/(2*dx))*(VORT[i+1][j]-VORT[i-1][j])
    return (J_1)

def J_2_func (PHI, VORT):
    J_2 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_2[i][j] = (1/(2*dx))*(PHI[i+1][j]*(1/(2*dy))*(VORT[i+1][j+1]-VORT[i+1][j-1]) - PHI[i-1][j]*(1/(2*dy))*(VORT[i-1][j+1]-VORT[i-1][j-1])) - (1/(2*dy))*(PHI[i][j+1]*(1/(2*dx))*(VORT[i+1][j+1]-VORT[i-1][j+1]) - PHI[i][j-1]*(1/(2*dx))*(VORT[i+1][j-1]-VORT[i-1][j-1]))    
    return (J_2)


def J_3_func (PHI, VORT):
    J_3 =  np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_3[i][j] = (1/(2*dy))*(VORT[i][j+1]*(1/(2*dx))*(PHI[i+1][j+1]-PHI[i-1][j+1]) - VORT[i][j-1]*(1/(2*dx))*(PHI[i+1][j-1] - PHI[i-1][j-1])) - (1/(2*dx))*(VORT[i+1][j]*(1/(2*dy))*(PHI[i+1][j+1]-PHI[i+1][j-1]) - VORT[i-1][j]*(1/(2*dy))*(PHI[i-1][j+1]-PHI[i-1][j-1]))
    return (J_3)

def J_final (PHI,VORT):
    return ((1/3)*(J_3_func(PHI,VORT)+J_2_func(PHI,VORT)+J_1_func(PHI,VORT)))



def MaxCourant(PHI, dt,dx ): #dx = 20
    absmaxU = np.amax(np.absolute(U_func(PHI)[2:-2,2:-2]))
    print("Abs max U is", absmaxU)
    absmaxV = np.amax(np.absolute(V_func(PHI)[2:-2,2:-2]))
    print("Abs max V is", absmaxV)
    return (absmaxU*dt/dx + absmaxV*dt/dx)

#IMPORTANT, RANGE (ghost nodes + walls ignored)
def totalKE (PHI):
    total = 0
    U = np.copy(U_func(PHI))
    V = np.copy(V_func(PHI))
    for i in range(2,Nx):
        for j in range(2,Ny):
            total += 0.5*(U[i][j]**2 + V[i][j]**2)
    return (total)

#IMPORTANT, VORTICITY AT WALL???????
def enstrophy (VORT):
    total = 0
    for i in range(2,101):
        for j in range(2,101):
            total += (VORT[i][j])**2
    return (total)
#NUMERICAL VORTICITY
def VORT_NUM_func (PHI):
    VORT_NUM = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1, Nx+1):
        for j in range(1,Ny+1):
            VORT_NUM[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j])+(1/dy**2)*(PHI[i][j-1]-2*PHI[i][j]+ PHI[i][j+1])
    return (VORT_NUM)
        

#******************************************************************************#
#***********************Initialize Vorticity***********************************#
#******************************************************************************#


VORT = np.zeros((Nx+2, Ny+2), dtype = float)

circ_factor = 8*10**(-5)

for i in range(1,Nx+1): #101 physical node
    for j in range(1,Ny+1):
        VORT[i][j] = circ_factor*math.exp(-(((i-1)*dx - x_int_1)**2 + ((j-1)*dy-y_int_1)**2)/a**2) + circ_factor*math.exp(-(((i-1)*dx - x_int_2)**2 + ((j-1)*dy-y_int_2)**2)/a**2)

        
#Boundary Condition

for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    VORT[i][1] = VORT[i][Ny]




print("finish VORT initialize")

#print(VORT[1:-1,1:-1])


plt.figure(0)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.title("Initial Vorticity")
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
#plt.savefig("INT,UNNORM-VORT-run_B.svg")






#******************************************************************************#
#*********************************** NORMALIZATION ****************************#
#******************************************************************************#

VORT_Avg = twoDSum(VORT[1:-1,1:-1])/(VORT[1:-1,1:-1].size)
VORT_loop_counter = 0

print("VORT loop counter", VORT_loop_counter)
print("mean vorticity = ", VORT_Avg)

"""
plt.figure(1)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
"""
#plt.savefig("INT,NORM-VORT-run_B.svg")
#VORT = np.transpose(VORT_R)


#******************************************************************************#
#*****************POISSON ADVANCING LISTS**************************************#
#******************************************************************************#

PHI = np.zeros((Nx+2, Ny+2), dtype = float)
Residual = np.zeros((Nx+2, Ny+2), dtype = float)


U = np.zeros((Nx+2, Ny+2), dtype = float)
V = np.zeros((Nx+2, Ny+2), dtype = float)
#Initial residual:             
for i in range( 0,Nx+2): 
    for j in range(0,Ny+2):
        Residual[i][j] =  -(VORT[i][j])

plt.figure(100)
plt.imshow( np.transpose(Residual[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()
plt.title("Initial Residual")
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")




"""
for a in range(0,Ny+2):
    Residual[Nx+1][a] = Residual[2][a]
    Residual[0][a] = Residual[Nx-1][a]
    Residual[1][a] = Residual[Nx][a]
for b in range(0,Nx+2):
    Residual[b][Ny+1] = Residual[b][2]
    Residual[b][0] = Residual[b][Ny-1]
    Residual[b][1] = Residual[b][Ny]

"""    
  
print("finish Residual initialization")


#test #print(u_cur.shape)
#result = 53x53
#*************************************************************************#
###########################Criterion for Convergence#######################
#*************************************************************************#
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
print("This is initial EPS = " ,EPS)


#******************************************************************************#
#************************************ALPHA*************************************#
#******************************************************************************#

#beta = dx/dy
#F-PLANE APPROXIMATION
sigma = (1/(1+(dx/dy)**2))*(math.cos(math.pi/Nx) + math.cos(math.pi/Ny)*(dx/dy)**2)
alpha = 2/(1+math.sqrt(1-sigma**2))

print("sigma =", sigma)    
print("alpha =", alpha)

timeStep = []
EPScounter = []

#******************************************************************************#
#***********************ADVANCING LIST (JACOBIANS)*****************************#
#******************************************************************************#
J_older = np.zeros((Nx+2, Ny+2), dtype = float) #n-2
J_old = np.zeros((Nx+2, Ny+2), dtype = float)   #n-1
J_cur = np.zeros((Nx+2, Ny+2), dtype = float)   #n
J_new = np.zeros((Nx+2, Ny+2), dtype = float)   #n+1

VORT_new = np.zeros((Nx+2, Ny+2), dtype = float)



#******************************************************************************#
#*********************** 0th Solver        ************************************#
#******************************************************************************#

EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0
maincounter = 0
while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    
    maincounter +=1
    #timeStep += [maincounter]
    #EPScounter +=[EPS]

print("0th Step: SOR Iteration Counter =", maincounter )



U = np.copy(U_func(PHI))
V = np.copy(V_func(PHI))


plt.figure(10000)
plt.title("Velocity Plot original")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(U[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()  

VORT_NUM = np.copy(VORT_NUM_func(PHI))
plt.figure(10020)
plt.title("numerical vort init")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none',extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()  


#plt.savefig("NumVort_RUN_B.svg", dpi = 300)

#******************************************************************************#
#*****************K total, Ensthropy total, max Courant************************#
#******************************************************************************#


print("0th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))
#print("0th Enstropy =", enstrophy(VORT))
      
#******************************************************************************#
#*****************************FIRST TIME STEP**********************************#
#******************************************************************************#

#Jacobian
J_older = np.copy(J_final(PHI,VORT))

#Step in time
for i in range(0,Nx+2):
    for j in range(0,Ny+2):
        VORT[i][j] = VORT[i][j] - J_older[i][j]*dt

#Boundary Condition
for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    #VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    #VORT[i][1] = VORT[i][Ny]

plt.figure(50)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
plt.title("I'm VORT")                     #URN
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()

    
#####POISSON SOLVER######

#ORDER is ******* IMPORTANT
PHI = np.zeros((Nx+2, Ny+2), dtype = float)
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0


maincounter = 0

while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    maincounter +=1
    
print("1th Step: SOR Iteration Counter =", maincounter )        
 
#UPDATING QUANTITIES
print("First Iteration Results:::::")
print("1th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))
#print("1st Enstropy =", enstrophy(VORT))
     



plt.figure(1493)

plt.title("VORT_NUM after 1st EF")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)

VORT_NUM = np.copy(VORT_NUM_func(PHI))  
#URN
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)




#******************************************************************************#
#*****************************SECOND TIME STEP*********************************#
#******************************************************************************#


#Jacobian
J_old = np.copy(J_final(PHI,VORT))

#Step in time
for i in range(0,Nx+2):
    for j in range(0,Ny+2):
        VORT[i][j] = VORT[i][j] - J_old[i][j]*dt

#Boundary Condition
for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    #VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    #VORT[i][1] = VORT[i][Ny]

plt.figure(5012)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
plt.title("VORT (2nd EF)")                     #URN
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()

    
#####POISSON SOLVER######

#ORDER is ******* IMPORTANT
PHI = np.zeros((Nx+2, Ny+2), dtype = float)
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0

maincounter = 0

while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j]
        for k in range(0,Ny+2):
            PHI[0][k] = PHI[Nx-1][k]
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
    #update finish
    for t in range(0,Nx+2):
        PHI[t][0] = PHI[t][Ny-1]
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    maincounter +=1
    
print("1th Step: SOR Iteration Counter =", maincounter )        
 
#UPDATING QUANTITIES
print("First Iteration Results:::::")
print("2th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))

#print("2nd Enstropy =", enstrophy(VORT))


#print("Enstropy =", enstrophy(VORT_NUM(PHI)))       



plt.figure(14931)

plt.title("VORT_NUM after 2st EF")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)

VORT_NUM = np.copy(VORT_NUM_func(PHI))  
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)

KEList = []
#******************************************************************************#
#****************************** MAIN TIME LOOP ********************************#
#******************************************************************************#

remove_counter=0


#while(remove_counter < 300):
for T in range(0,Nt):
    
    #Jacobian
    J_cur = np.copy(J_final(PHI,VORT))

    #Step in time
    for i in range(0,Nx+2):
        for j in range(0,Ny+2):
            VORT[i][j] = VORT[i][j] - dt*(23/12)*J_cur[i][j] + dt*(16/12)*J_old[i][j] - dt*(5/12)*J_older[i][j]

    #Boundary Condition
    for j in range(0,Ny+2):
        VORT[Nx+1][j] = 0
        VORT[0][j] = 0
        #VORT[1][j] = VORT[Nx][j]
    for i in range(0,Nx+2):
        VORT[i][Ny+1] = 0
        VORT[i][0] = 0






    #####POISSON SOLVER######

    #ORDER is ******* IMPORTANT
    PHI = np.zeros((Nx+2, Ny+2), dtype = float)
    EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    # = 1.0

    maincounter = 0

    while (EPS > tol):
        for j in range(2,Ny):  
            #update ghost nodes along phi[i][Ny+1]
            if (j == Ny):  #IMPORTANT , NY-1 for wall
                for z in range(0,Nx+2): 
                    PHI[z][Ny+1] = PHI[z][2]
            #update finish
            for i in range(2,Nx): #101 physical node WALL
                #update ghost nodes along phi[Nx+1][j]
                if (i == Nx):
                    for k in range(0,Ny+2):
                        PHI[Nx+1][k] = PHI[2][k]
                #update finisH
                Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
                PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))

            #Update ghosts at PHI[0][j]
            for k in range(0,Ny+2):
                PHI[0][k] = PHI[Nx-1][k]
                PHI[1][k] = PHI[Nx][k] #make sure left right equal

        #update finish
        for t in range(0,Nx+2):
            PHI[t][0] = PHI[t][Ny-1]
            PHI[t][1] = PHI[t][Ny] 
        #make sure left right equal

        #update finish
        #Update ghosts at PHI[i][0]
        tempResidual = np.copy(Residual)
        tempResidual = np.absolute(tempResidual[2:-2,2:-2])

        EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
        #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
        #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))

        maincounter +=1

    print("big counter =", remove_counter)
    remove_counter+=1
    
    #Shifting:
    J_older = np.copy(J_old)
    J_old = np.copy(J_cur)
    

    if (remove_counter%10 == 0):
        
        plt.figure(remove_counter+102432)
        x = np.linspace(0,2000,101)
        y = np.linspace(0,2000,101)
        plt.xlabel("X Domain [m]")
        plt.ylabel("Y Domain [m]")
        X,Y = np.meshgrid(x,y)
        plt.title("VORT")                     #URN
        plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
        plt.colorbar()
        plt.savefig("VORT"+str(remove_counter)+".svg")
        


        plt.figure(remove_counter+102433)
        x = np.linspace(0,2000,101)
        y = np.linspace(0,2000,101)
        plt.xlabel("X Domain [m]")
        plt.ylabel("Y Domain [m]")
        X,Y = np.meshgrid(x,y)
        plt.title("PHI")                     #URN
        plt.imshow( np.transpose(PHI[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
        plt.colorbar()
        plt.savefig("PHI"+str(remove_counter)+".svg")
        print("remove_counter =", remove_counter)
        print("nth Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
        print("Total KE =", totalKE(PHI))
        


plt.figure(475)
time = np.linspace(0,1728000,len(KEList))
plt.plot(time, KEList)


# In[ ]:





# In[1]:


#*****************NOTES:********************#

"""
Most update version (SAT 12:00am)
RUN: Non-inverted Original Run



Max Courant printout silenced
Have list to store KE, MaxCourtant & Ensthropy




"""


# from numpy import math
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import matplotlib
import math
import pandas as pd
np.core
from sys import maxsize
from numpy import set_printoptions
set_printoptions(threshold=maxsize)

print("hi")

#*****************Constants********************#
Lx = 2000*1000
Nx = 101
Ly = 2000*1000
Ny = 101
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
print("dx = ",dx)
print("dy = ",dy)


#r = 500.0
tol = 10**-8

a = 180*1000
b = 600*1000

x_int_1 = (Lx/2) + (b/2)
y_int_1 = Ly/2
x_int_2 = (Lx/2) - (b/2)
y_int_2 = Ly/2

Ttot = 1728000 #SECONDS
dt = 500
Nt = int(Ttot/dt)


beta = 0

#****************************************#
#*************Functions******************#
#****************************************#

def forceEqual(array): #left right symmetry
    finarr = np.transpose(array)
    i = 0
    while (i != 52):
        array[i] = array[102-i]
        i+=1
    return np.transpose(finarr)

def twoDSum (array):
    newarr = []
    for i in array:
        newarr+=[math.fsum(i)]
    return (math.fsum(newarr))


def U_func (PHI):
    dy = 20000
    U = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            U[i][j] = (PHI[i][j-1]-PHI[i][j+1])/(2*dy)
    return (U)
    
def V_func (PHI):
    dx = 20000
    V = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            V[i][j] = (PHI[i+1][j]-PHI[i-1][j])/(2*dx)
    return (V)


def J_1_func (PHI, VORT): #return standard size Nx+2 x Ny+2 array
    J_1 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1, Nx+1):
        for j in range(1,Ny+1):
            J_1[i][j] = (1/(2*dx))*(PHI[i+1][j]-PHI[i-1][j])*(1/(2*dy))*(VORT[i][j+1] - VORT[i][j-1]) - (1/(2*dy))*(PHI[i][j+1] - PHI[i][j-1])*(1/(2*dx))*(VORT[i+1][j]-VORT[i-1][j])
    return (J_1)

def J_2_func (PHI, VORT):
    J_2 = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_2[i][j] = (1/(2*dx))*(PHI[i+1][j]*(1/(2*dy))*(VORT[i+1][j+1]-VORT[i+1][j-1]) - PHI[i-1][j]*(1/(2*dy))*(VORT[i-1][j+1]-VORT[i-1][j-1])) - (1/(2*dy))*(PHI[i][j+1]*(1/(2*dx))*(VORT[i+1][j+1]-VORT[i-1][j+1]) - PHI[i][j-1]*(1/(2*dx))*(VORT[i+1][j-1]-VORT[i-1][j-1]))    
    return (J_2)


def J_3_func (PHI, VORT):
    J_3 =  np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1,Nx+1):
        for j in range(1,Ny+1):
            J_3[i][j] = (1/(2*dy))*(VORT[i][j+1]*(1/(2*dx))*(PHI[i+1][j+1]-PHI[i-1][j+1]) - VORT[i][j-1]*(1/(2*dx))*(PHI[i+1][j-1] - PHI[i-1][j-1])) - (1/(2*dx))*(VORT[i+1][j]*(1/(2*dy))*(PHI[i+1][j+1]-PHI[i+1][j-1]) - VORT[i-1][j]*(1/(2*dy))*(PHI[i-1][j+1]-PHI[i-1][j-1]))
    return (J_3)

def J_final (PHI,VORT):
    return ((1/3)*(J_3_func(PHI,VORT)+J_2_func(PHI,VORT)+J_1_func(PHI,VORT)))



def MaxCourant(PHI, dt,dx ): #dx = 20
    absmaxU = np.amax(np.absolute(U_func(PHI)[2:-2,2:-2]))
    #print("Abs max U is", absmaxU)
    absmaxV = np.amax(np.absolute(V_func(PHI)[2:-2,2:-2]))
    #print("Abs max V is", absmaxV)
    return (absmaxU*dt/dx + absmaxV*dt/dx)

#IMPORTANT, RANGE (ghost nodes + walls ignored)
def totalKE (PHI):
    total = 0
    U = np.copy(U_func(PHI))
    V = np.copy(V_func(PHI))
    for i in range(2,Nx):
        for j in range(2,Ny):
            total += 0.5*(U[i][j]**2 + V[i][j]**2)
    return (total)

#IMPORTANT, VORTICITY AT WALL???????
def enstrophy (VORT):
    total = 0
    for i in range(2,101):
        for j in range(2,101):
            total += (VORT[i][j])**2
    return (total)

#NUMERICAL VORTICITY
def VORT_NUM_func (PHI):
    VORT_NUM = np.zeros((Nx+2, Ny+2), dtype = float)
    for i in range(1, Nx+1):
        for j in range(1,Ny+1):
            VORT_NUM[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j])+(1/dy**2)*(PHI[i][j-1]-2*PHI[i][j]+ PHI[i][j+1])
    return (VORT_NUM)
        

#******************************************************************************#
#***********************Initialize Vorticity***********************************#
#******************************************************************************#


VORT = np.zeros((Nx+2, Ny+2), dtype = float)

circ_factor = 8*10**(-5)

for i in range(1,Nx+1): #101 physical node
    for j in range(1,Ny+1):
        VORT[i][j] = circ_factor*math.exp(-(((i-1)*dx - x_int_1)**2 + ((j-1)*dy-y_int_1)**2)/a**2) + circ_factor*math.exp(-(((i-1)*dx - x_int_2)**2 + ((j-1)*dy-y_int_2)**2)/a**2)

        
#Boundary Condition

for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    VORT[i][1] = VORT[i][Ny]




print("finish VORT initialize")

#print(VORT[1:-1,1:-1])


plt.figure(0)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.title("Initial Vorticity")
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
#plt.savefig("INT,UNNORM-VORT-run_B.svg")






#******************************************************************************#
#*********************************** NORMALIZATION ****************************#
#******************************************************************************#

VORT_Avg = twoDSum(VORT[1:-1,1:-1])/(VORT[1:-1,1:-1].size)
VORT_loop_counter = 0

print("VORT loop counter", VORT_loop_counter)
print("mean vorticity = ", VORT_Avg)

"""
plt.figure(1)
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()

plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
"""
#plt.savefig("INT,NORM-VORT-run_B.svg")
#VORT = np.transpose(VORT_R)


#******************************************************************************#
#*****************POISSON ADVANCING LISTS**************************************#
#******************************************************************************#

PHI = np.zeros((Nx+2, Ny+2), dtype = float)
Residual = np.zeros((Nx+2, Ny+2), dtype = float)



U = np.zeros((Nx+2, Ny+2), dtype = float)
V = np.zeros((Nx+2, Ny+2), dtype = float)
#Initial residual:             
for i in range( 0,Nx+2): 
    for j in range(0,Ny+2):
        Residual[i][j] =  -(VORT[i][j])

plt.figure(100)
plt.imshow( np.transpose(Residual[1:-1,1:-1]), interpolation='none', extent=[0,100,0,100])
plt.colorbar()
plt.title("Initial Residual")
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")




"""
for a in range(0,Ny+2):
    Residual[Nx+1][a] = Residual[2][a]
    Residual[0][a] = Residual[Nx-1][a]
    Residual[1][a] = Residual[Nx][a]
for b in range(0,Nx+2):
    Residual[b][Ny+1] = Residual[b][2]
    Residual[b][0] = Residual[b][Ny-1]
    Residual[b][1] = Residual[b][Ny]

"""    
  
print("finish Residual initialization")


#test #print(u_cur.shape)
#result = 53x53
#*************************************************************************#
###########################Criterion for Convergence#######################
#*************************************************************************#
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
print("This is initial EPS = " ,EPS)


#******************************************************************************#
#************************************ALPHA*************************************#
#******************************************************************************#


#F-PLANE APPROXIMATION


sigma = (1/(1+(dx/dy)**2))*(math.cos(math.pi/Nx) + math.cos(math.pi/Ny)*(dx/dy)**2)
alpha = 2/(1+math.sqrt(1-sigma**2))

print("sigma =", sigma)    
print("alpha =", alpha)

timeStep = []
EPScounter = []

#******************************************************************************#
#***********************ADVANCING LIST (JACOBIANS)*****************************#
#******************************************************************************#
J_older = np.zeros((Nx+2, Ny+2), dtype = float) #n-2
J_old = np.zeros((Nx+2, Ny+2), dtype = float)   #n-1
J_cur = np.zeros((Nx+2, Ny+2), dtype = float)   #n
J_new = np.zeros((Nx+2, Ny+2), dtype = float)   #n+1


PHI_older = np.zeros((Nx+2, Ny+2), dtype = float) #n-2
PHI_old= np.zeros((Nx+2, Ny+2), dtype = float)   #n-1
PHI_cur = np.zeros((Nx+2, Ny+2), dtype = float)#n
PHI_new = np.zeros((Nx+2, Ny+2), dtype = float)
#VORT_new = np.zeros((Nx+2, Ny+2), dtype = float)


KE_List = []
Ens_List = []
Cour_List = []
#******************************************************************************#
#*********************** 0th Solver        ************************************#
#******************************************************************************#


for i in range(2,Nx):
    for j in range(2,Ny):
        Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
        
EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0
maincounter = 0
while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j] NEWWW
        for k in range(0,Ny+2):
            PHI[1][k] = 0
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
 
    #update finish
    for t in range(0,Nx+2):
        PHI[t][Nx+1] = 0
        PHI[t][0] = 0
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    
    maincounter +=1
    #timeStep += [maincounter]
    #EPScounter +=[EPS]

print("0th Step: SOR Iteration Counter =", maincounter )



U = np.copy(U_func(PHI))
V = np.copy(V_func(PHI))


plt.figure(10000)
plt.title("Velocity Plot original")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(U[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()  

VORT_NUM = np.copy(VORT_NUM_func(PHI))
plt.figure(10020)
plt.title("numerical vort init")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
X,Y = np.meshgrid(x,y)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none',extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()  


#plt.savefig("NumVort_RUN_B.svg", dpi = 300)

#******************************************************************************#
#*****************K total, Ensthropy total, max Courant************************#
#******************************************************************************#


print("0th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))
print("0th Enstropy =", enstrophy(VORT))
KE_List += [totalKE(PHI)]
Ens_List+= [enstrophy(VORT)]
Cour_List += [MaxCourant(PHI, dt, dx )]    

#******************************************************************************#
#*****************************FIRST TIME STEP**********************************#
#******************************************************************************#

#Jacobian
J_older = np.copy(J_final(PHI,VORT))

PHI_older = np.copy(PHI)

#Step in time
for i in range(2,Nx):
    for j in range(2,Ny):
        VORT[i][j] = VORT[i][j] - J_older[i][j]*dt - beta*(PHI_older[i+1][j]-PHI_older[i-1][j])/(2*dx)

#Boundary Condition
for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    #VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    #VORT[i][1] = VORT[i][Ny]

plt.figure(50)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
plt.title("VORT")                     #URN
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()

    
#####POISSON SOLVER######

#ORDER is ******* IMPORTANT
#PHI = np.zeros((Nx+2, Ny+2), dtype = float)

for i in range(2,Nx):
    for j in range(2,Ny):
        Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]

EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0


maincounter = 0


while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j] NEWWW
        for k in range(0,Ny+2):
            PHI[1][k] = 0
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
 
    #update finish
    for t in range(0,Nx+2):
        PHI[t][Nx+1] = 0
        PHI[t][0] = 0
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    
    maincounter +=1
    
print("1th Step: SOR Iteration Counter =", maincounter )        
 
#UPDATING QUANTITIES
print("First Iteration Results:::::")
print("1th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))
print("1st Enstropy =", enstrophy(VORT))
KE_List += [totalKE(PHI)]
Ens_List+= [enstrophy(VORT)]
Cour_List += [MaxCourant(PHI, dt, dx )]  
     



plt.figure(1493)

plt.title("VORT_NUM after 1st EF")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)

VORT_NUM = np.copy(VORT_NUM_func(PHI))  
#URN
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)




#******************************************************************************#
#*****************************SECOND TIME STEP*********************************#
#******************************************************************************#


#Jacobian
J_old = np.copy(J_final(PHI,VORT))
PHI_old = np.copy(PHI)

#Step in time
for i in range(2,Nx):
    for j in range(2,Ny):
        VORT[i][j] = VORT[i][j] - J_old[i][j]*dt - beta*(PHI_old[i+1][j]-PHI_old[i-1][j])/(2*dx)

#Boundary Condition
for j in range(0,Ny+2):
    VORT[Nx+1][j] = 0
    VORT[0][j] = 0
    #VORT[1][j] = VORT[Nx][j]
for i in range(0,Nx+2):
    VORT[i][Ny+1] = 0
    VORT[i][0] = 0
    #VORT[i][1] = VORT[i][Ny]

plt.figure(5012)
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)
plt.title("VORT (2nd EF)")                     #URN
plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()

    
#####POISSON SOLVER######

#ORDER is ******* IMPORTANT
#PHI = np.zeros((Nx+2, Ny+2), dtype = float)


for i in range(2,Nx):
    for j in range(2,Ny):
        Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]

EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
# = 1.0

maincounter = 0

while (EPS > tol):
    for j in range(2,Ny):  
        #update ghost nodes along phi[i][Ny+1]
        if (j == Ny):  #IMPORTANT , NY-1 for wall
            for z in range(0,Nx+2): 
                PHI[z][Ny+1] = PHI[z][2]
        #update finish
        for i in range(2,Nx): #101 physical node WALL
            #update ghost nodes along phi[Nx+1][j]
            if (i == Nx):
                for k in range(0,Ny+2):
                    PHI[Nx+1][k] = PHI[2][k]
            #update finisH
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
            PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))
        
        #Update ghosts at PHI[0][j] NEWWW
        for k in range(0,Ny+2):
            PHI[1][k] = 0
            PHI[1][k] = PHI[Nx][k] #make sure left right equal
            
 
    #update finish
    for t in range(0,Nx+2):
        PHI[t][Nx+1] = 0
        PHI[t][0] = 0
        PHI[t][1] = PHI[t][Ny] 
    #make sure left right equal
    
    #update finish
    #Update ghosts at PHI[i][0]
    tempResidual = np.copy(Residual)
    tempResidual = np.absolute(tempResidual[2:-2,2:-2])
    
    EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
    #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))
    
    
    maincounter +=1
    
print("1th Step: SOR Iteration Counter =", maincounter )        
 
#UPDATING QUANTITIES
print("First Iteration Results:::::")
print("2th Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
print("Total KE =", totalKE(PHI))
print("2nd Enstropy =", enstrophy(VORT))

KE_List += [totalKE(PHI)]
Ens_List+= [enstrophy(VORT)]
Cour_List += [MaxCourant(PHI, dt, dx )]  


#print("Enstropy =", enstrophy(VORT_NUM(PHI)))       



plt.figure(14931)

plt.title(" VORT_NUM after 2st EF")
x = np.linspace(0,2000,101)
y = np.linspace(0,2000,101)
plt.xlabel("X Domain [m]")
plt.ylabel("Y Domain [m]")
X,Y = np.meshgrid(x,y)

VORT_NUM = np.copy(VORT_NUM_func(PHI))  
plt.imshow( np.transpose(VORT_NUM[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
plt.colorbar()
#plt.savefig("NumVort_RUN_B.svg", dpi = 300)


#******************************************************************************#
#****************************** MAIN TIME LOOP ********************************#
#******************************************************************************#

remove_counter=0


#while(remove_counter < 300):
for T in range(0,Nt):
    
    #Jacobian
    J_cur = np.copy(J_final(PHI,VORT))
    PHI_cur = np.copy(PHI)
    
    #Step in time
    for i in range(2,Nx):
        for j in range(2,Ny):
            VORT[i][j] = VORT[i][j] - (23/12)*(dt*J_cur[i][j] + beta*(PHI_cur[i+1][j]-PHI_cur[i-1][j])/(2*dx)) + (16/12)*(J_old[i][j]*dt + beta*(PHI_old[i+1][j]-PHI_old[i-1][j])/(2*dx)) - (5/12)*(dt*J_older[i][j] + beta*(PHI_older[i+1][j]-PHI_older[i-1][j])/(2*dx))

    #Boundary Condition
    for j in range(0,Ny+2):
        VORT[Nx+1][j] = 0
        VORT[0][j] = 0
        #VORT[1][j] = VORT[Nx][j]
    for i in range(0,Nx+2):
        VORT[i][Ny+1] = 0
        VORT[i][0] = 0






    #####POISSON SOLVER######

    #ORDER is ******* IMPORTANT
    #PHI = np.zeros((Nx+2, Ny+2), dtype = float)
    
    
    for i in range(2,Nx):
        for j in range(2,Ny):
            Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]

              
    
    EPS = np.amax(np.absolute(Residual))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
    # = 1.0
   
    maincounter = 0

    
    while (EPS > tol):
        for j in range(2,Ny):  
            #update ghost nodes along phi[i][Ny+1]
            if (j == Ny):  #IMPORTANT , NY-1 for wall
                for z in range(0,Nx+2): 
                    PHI[z][Ny+1] = PHI[z][2]
            #update finish
            for i in range(2,Nx): #101 physical node WALL
                #update ghost nodes along phi[Nx+1][j]
                if (i == Nx):
                    for k in range(0,Ny+2):
                        PHI[Nx+1][k] = PHI[2][k]
                #update finisH
                Residual[i][j] = (1/dx**2)*(PHI[i-1][j] -2*PHI[i][j] + PHI[i+1][j]) + (1/dy**2)*(PHI[i][j-1] -2*PHI[i][j] + PHI[i][j+1]) - VORT[i][j]
                PHI[i][j] = PHI[i][j] + (alpha*Residual[i][j])/(2*((1/dx**2)+(1/dy**2)))

            #Update ghosts at PHI[0][j] NEWWW
            for k in range(0,Ny+2):
                PHI[1][k] = 0
                PHI[Nx][k] = 0 #make sure left right equal


        #update finish
        for t in range(0,Nx+2):
            PHI[t][1] = 0
            PHI[t][Ny] = 0
        #make sure left right equal

        #update finish
        #Update ghosts at PHI[i][0]
        tempResidual = np.copy(Residual)
        tempResidual = np.absolute(tempResidual[2:-2,2:-2])

        EPS = np.amax(np.absolute(Residual[2:-2,2:-2]))/((2*((1/dx**2)+(1/dy**2))*twoDSum(np.absolute(PHI[2:-2,2:-2])))+ np.amax(np.absolute(VORT[2:-2,2:-2])))
        #print('EPS is', EPS, '  MainCounter =', maincounter, "Max of abs Value = "  , np.amax(np.absolute(Residual[2:-2,2:-2])))
        #print(" Index of max abs = ", np.unravel_index(tempResidual.argmax(), tempResidual.shape))


        maincounter +=1

    print("big counter =", remove_counter, "converge counter =", maincounter)
    remove_counter+=1
    
    #output important quantities
    KE_List += [totalKE(PHI)]
    Ens_List+= [enstrophy(VORT)]
    Cour_List += [MaxCourant(PHI, dt, dx )] 
    
    
    #Shifting:
    J_older = np.copy(J_old)
    J_old = np.copy(J_cur)
    PHI_older = np.copy(PHI_old)
    PHI_old = np.copy(PHI_cur)
    
    

    if (remove_counter%25 == 0):
        
        plt.figure(remove_counter+102432)
        x = np.linspace(0,2000,101)
        y = np.linspace(0,2000,101)
        plt.xlabel("X Domain [m]")
        plt.ylabel("Y Domain [m]")
        X,Y = np.meshgrid(x,y)
        plt.title("VORT")                     #URN
        plt.imshow( np.transpose(VORT[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
        plt.colorbar()



        plt.figure(remove_counter+102433)
        x = np.linspace(0,2000,101)
        y = np.linspace(0,2000,101)
        plt.xlabel("X Domain [m]")
        plt.ylabel("Y Domain [m]")
        X,Y = np.meshgrid(x,y)
        plt.title(" PSI ")                     #URN
        plt.imshow( np.transpose(PHI[1:-1,1:-1]), interpolation='none', extent=[0,Lx/1000,0,Ly/1000])
        plt.colorbar()
        print("remove_counter =", remove_counter)
        print("nth Step: Max Courant = ", MaxCourant(PHI, dt, dx ))
        print("Total KE =", totalKE(PHI))
        print("Ensthropy", enstrophy(VORT))

        
        #visualize difference (plot difference between both vorticity)
        
        
        
        
KE_List += [totalKE(PHI)]
Ens_List+= [enstrophy(VORT)]
Cour_List += [MaxCourant(PHI, dt, dx )]          
    
plt.figure(475)
plt.title("Total Kinetic Energy")
x = np.linspace(1,len(KE_List),len(KE_List))
plt.plot(x,KE_List)


plt.figure(476)
plt.title("Enstrophy")
x = np.linspace(1,len(Ens_List),len(Ens_List))
plt.plot(x,Ens_List)

plt.figure(477)
plt.title("Maximum Courant Number")
x = np.linspace(1,len(Cour_List),len(Cour_List))
plt.plot(x,Cour_List)



# In[6]:


plt.figure(470)
#plt.title("Maximum Courant Number")
#plt.xscale("log")
#plt.yscale("log")
x = np.linspace(1,len(Ens_List),len(Ens_List))
plt.plot(x,Ens_List)
plt.ylabel("Enstrophy")
plt.xlabel("Time Step")


# In[ ]:




