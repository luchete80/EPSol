# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
from numpy import *
import numpy.matlib

import array as arr
#Data
lx=1.
ly=1.
nex=1
ney=1
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

#X is combines per row to calculate J
p=1.0/1.732050807568877
gauss=[-p,p]
wi=[0.777,0.777]
#Numerated as in Bathe
X2=numpy.matlib.zeros((4, 2))
#Numerated as in deal.ii
#X2=matrix([[0,0],[1,0],[0,1],[1,1]])

K=matrix(numpy.matlib.zeros((2*numnodes, 2*numnodes)))
R=matrix(numpy.matlib.zeros((2*numnodes, 1)))

Kel=matrix(numpy.matlib.zeros((8, 8)))
B  =matrix(numpy.matlib.zeros((3, 8)))
dHxy=matrix(numpy.matlib.zeros((2, 4)))

ey=200.e9
nu=0.33

ck = ey*(1. - nu) / ((1. + nu)*(1. - 2. * nu))
c=matrix(numpy.matlib.zeros((3, 3)))
c[0,0]=c[1,1]=ck;					
c[0,1]=c[1,0]=ck*nu / (1. - nu)
c[2,2]=ck*(1 - 2 * nu) / (2.*(1. - nu))
print ("C matrix")
print(c)
for e in range (numel):
    Kel=0.
    for n in range(4):
        X2[n]=node[elnodes.astype(int)[e][n]]
    print ("Element Nodes")
    print (X2)
    for i in range(2):
        for j in range(2):
            r=gauss[i]
            s=gauss[j]

            #Numerated as in Bathe
            dHrs=matrix([[(1+s),-(1+s),-(1-s),(1-s)], [(1+r),(1-r),-(1-r),-(1+r)] ])
            #Numerated as in deal.ii
            #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
            dHrs/=4
            J=dHrs*X2
            dHxy=linalg.inv(J)*dHrs
            for k in range(4):
                B[0,2*k  ]=dHxy[0,k]
                B[1,2*k+1]=dHxy[1,k]
                B[2,2*k  ]=dHxy[1,k]
                B[2,2*k+1]=dHxy[0,k]
            w=linalg.det(J)*1.
            # print ("weight")
            # print (w)
            #print (B)
            Kel+=(B.transpose()*c*B*w)
            #print (K)
    print(Kel)

    #Assembly Matrix
    for n in range (4):
        d=elnodes.astype(int)[e][n]
        vn[2*n  ]=2*d
        vn[2*n+1]=2*d+1
    print(vn.astype(int))
    
    for row in range(8):
        for col in range (8):
            K[vn.astype(int)[row],vn.astype(int)[col]]=Kel[row,col]
#print (K)

#Boundary conditions
#Bottom left node both directions
#Bottom right node, vertical direction
fix=[0,1,2*(nex+1)-1]
for f in range (3):
    dof=fix[f]
    for i in range (2*numnodes):
        K[dof,i]=0.
        K[i,dof]=0.
        K[dof,dof]=1.

#Upper left node, x direction force			
R[(2*nex)*(2*nex),0]=1000.
print(R)

U=linalg.solve(K, R)
print ("Results")
print(U)
