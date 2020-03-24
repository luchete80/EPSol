# https://docs.sympy.org/latest/tutorial/matrices.html
# From maniatty rigid visco plastic
# Stabilized finite element method for viscoplastic ow:
# formulation with state variable evolution

# Plain strain/stress 1 element example
#import math
from numpy import *
import numpy.matlib #for zeros

import array as arr

#import deriv

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

X2=numpy.matlib.zeros((4, 2))
p=1.0/1.732050807568877
gauss=[-p,p]

#DOFS
var_dim =[2,1,1]
var_edof=[8,4,4]
var_edof = array(var_edof, dtype=int)

#Shape functions
Nps=matrix(numpy.matlib.zeros((1, 4)))
Nv=matrix(numpy.matlib.zeros((2, 8)))

#Matrices
BvtBv=matrix(numpy.matlib.zeros((nodxel*dim, nodxel*dim)))
temp8x8=matrix(numpy.matlib.zeros((nodxel*dim, nodxel*dim)))
Bv=matrix(numpy.matlib.zeros((4, nodxel*dim)))
Bps=matrix(numpy.matlib.zeros((2, 8)))

UV  =matrix(numpy.matlib.zeros((8, 1)))
Up  =matrix(numpy.matlib.zeros((4, 1)))
Us  =matrix(numpy.matlib.zeros((4, 1)))

Kt=[ [matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[1])))  
     ]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[1]))),
     ]
    ]     

#Kt
# Elemental matrix
# 16x16
# (8x8) (8x4) (8x4) 
# (4x8) (4x4) (4x4)
# (4x8) (4x4) (4x4) 


end=0

while (end==0):
    end=1
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
                Nps  =0.25*matrix([(1+sg)*(1+rg),(1-rg)*(1+sg),(1-sg)*(1-rg),(1-sg)*(1+rg)])   
                dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
                #Numerated as in deal.ii
                #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
                dHrs/=4.
                J=dHrs*X2
                dHxy=linalg.inv(J)*dHrs
                detJ=linalg.det(J)
                #Calculate shape functions
                #Bs=J-1 dHrs(B.13)
                Bps=dHxy
                for k in range(nodxel):
                    #shape functions
                    Nv[0,2*k  ]=Nv[1,2*k+1]=Nps[0,k]
                    
                dp=Nps*Up
                ds=Nps*Us
                tau=1.
                he=1.
                alpha=1.
                
                #p+=dp
                #s+=ds
                w=1.
                wJ=w*detJ                
                #

                BvtBv=Bv.transpose*Bv
                for i in range(8):
                    temp8x8[i,i]=BbtBv
                Kt[0][0]=Kt[0][0]+(
                            temp8x8+BvtBv
                            )*wJ
                
                #Kt u p (8x4) div * N p 
                #b(dp,v*)=-dp x div(v*)
                Kt[0][1]=Kt[0][1]-Bv.transpose()*Nps
                print ("Kt 0 1 ",Kt[0][1])
     
                print ("Kt 0 2 ",Kt[0][2])                
                Kt[1][0]=Kt[1][0]+(
                            *Nps.transpose() #Div dv 
                            )*wJ
                
                #mstab=tau[grad(dp).grad(p*)] 
                print("Kt11",Kt[1][1])
                #K11 is (4x4)
                Kt[1][1]=   Kt[1][1]+(
                            float(alpha*he**2.)*Bps.transpose()*Bps
                            )*wJ
