# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
#import math
from numpy import *
import numpy.matlib #for zeros

import array as arr

import deriv

#----------------------------
#Input Data------------------
form=2
lx=1.
ly=1.
nex=2
ney=2

nodxel=4
dim=2
#-------------


#***** MESH**********
dx=lx/nex
dy=ly/ney
numel=nex*ney
vn=zeros(8) #Dof connectivity
numnodes=(nex+1)*(ney+1)
node=zeros((numnodes, 2))
elnodes=zeros((numel, 4))#Connectivity
elnodes.astype(int)
#print (node)
#Mesh generation 
y=0.
n=0
for nx in range (ney+1):
    x=0.
    for ny in range (nex+1):
        #print("Node: "+str(x)+" , "+str(y))
        #this->node.push_back(Node(n,x,y,0.0))
        node[n]=[x,y]
        x=x+dx
        n=n+1
    y=y+dy
    
print(node)
#Connectivity
e=0
for ey in range (ney):
    for ex in range (nex):
        elnodes[e]=[(nex+1)*(ey+1)+ex+1,
                    (nex+1)*(ey+1)+ex,
                    (nex+1)*ey + ex,
                    (nex+1)*ey + ex+1]
        e=e+1
print(elnodes)
#-------------------------- MESH
#Formulation type and DOFs
if form==1:
    ndof =11
else:
    ndof =12
    
ey=200.e9
nu=0.33

ck = ey/ ((1. + nu)*(1. - 2. * nu))
c=matrix(numpy.matlib.zeros((4, 4)))
c[0,0]=c[1,1]=c[2,2]=ck*(1.-nu)					
c[0,1]=c[1,0]=c[0,2]=c[2,0]=c[2,1]=c[1,2]=ck*nu
c[3,3]=ck*(.5 - nu)

#Matrices
Bv=matrix(numpy.matlib.zeros((4, nodxel*dim)))
Nv=matrix(numpy.matlib.zeros((dim, nodxel*dim)))

for e in range (4):
    #Obtain Ve from global
    Kel=0.
    for n in range(4):
        X2[n]=node[elnodes.astype(int)[e][n]]
    print ("Element ", e)
    print ("Element Nodes")
    print (X2)
    
    print ("Nv",Nv)
    for ig in range(2):
        for jg in range(2):
            rg=gauss[ig]
            sg=gauss[jg]

            #Numerated as in Bathe
            Ns  =0.25*matrix([(1+sg)*(1+rg),(1-rg)*(1+sg),(1-sg)*(1-rg),(1-sg)*(1+rg)])   
            dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
            #Numerated as in deal.ii
            #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
            dHrs/=4.
            J=dHrs*X2
            dHxy=linalg.inv(J)*dHrs
            detJ=linalg.det(J)
            #Calculate shape functions
            #Bs=J-1 dHrs(B.13)
            Bs=dHxy
            for k in range(nodxel):
                #shape functions
                Nv[0,2*k  ]=Nv[1,2*k+1]=Ns[0,k]
                
 
            #From Zienkiewicz page 285 11.11
            Ka=Bv.transpose()*Dd*B*wJ
            Kc=
            Kp=Np.transpose()*float(1.0/k_pen)*Np*wJ