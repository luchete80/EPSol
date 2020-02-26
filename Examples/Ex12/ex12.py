# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
from numpy import *
import numpy.matlib #for zeros

import array as arr

#X is combines per row to calculate J
p=1.0/1.732050807568877
gauss=[-p,p]
wi=[0.777,0.777]
#Numerated as in Bathe
X2=matrix([[1,1],[0,1],[0,0],[1,0]])
#Numerated as in deal.ii
#X2=matrix([[0,0],[1,0],[0,1],[1,1]])

#Nodal values
Ve=matrix(numpy.matlib.zeros((1, 8)))

#Shape Functions
Ns=matrix(numpy.matlib.zeros((1, 4)))
Nv=matrix(numpy.matlib.zeros((2, 8)))
NsigF=matrix(numpy.matlib.zeros((4, 16)))

#Derivatives
dHxy=matrix(numpy.matlib.zeros((2, 4)))
Bs=matrix(numpy.matlib.zeros((2, 8)))
Bv=matrix(numpy.matlib.zeros((4, 8)))
LM=matrix(numpy.matlib.zeros((4, 4)))

K=matrix(numpy.matlib.zeros((44, 44)))
R=matrix(numpy.matlib.zeros((44, 1)))
dVxy=matrix(numpy.matlib.zeros((4, 2)))


#Material properties (Table 2.1 p28)
#HSLA-65 steel with strain rate between e-3 and e-4
A0=6.34e11
psi=3.25
m=0.1956
n=0.06869
s0=80.0e6


print (K)
print (X2[0,0])
print (X2[2,1])

ey=200.e9
nu=0.33

ck = ey*(1. - nu) / ((1. + nu)*(1. - 2. * nu))
c=matrix(numpy.matlib.zeros((3, 3)))
c[0,0]=c[1,1]=ck;					
c[0,1]=c[1,0]=ck*nu / (1. - nu)
c[2,2]=ck*(1 - 2 * nu) / (2.*(1. - nu))
print ("C matrix")
print(c)
#Set G and H matrices
#2.32
G=matrix([[1,0,0,0],[0,0,0,1],[0,0,0,1],[0,1,0,0]]) 
H=matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]) 
    
for e in range (4):
    #Obtain Ve from global
    for i in range(2):
        for j in range(2):
            r=gauss[i]
            s=gauss[j]

            #Numerated as in Bathe
            Ns  =0.25*matrix([(1+s)*(1+r),(1-r)*(1+s),(1-s)*(1-r),(1-s)*(1-r)])            
            dHrs=matrix([[(1+s),-(1+s),-(1-s),(1-s)], [(1+r),(1-r),-(1-r),-(1+r)] ])
            #Numerated as in deal.ii
            #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
            dHrs/=4
            J=dHrs*X2
            dHxy=linalg.inv(J)*dHrs
            #Calculate shape functions
            #Bs=J-1 dHrs(B.13)
            Bs=dHxy
            for k in range(4):
                #shape functions
                Nv[0,2*k  ]=Nv[1,2*k+1]=Ns[0,k]
                #derivatives Bv (B.14)
                Bv[0,2*k  ]=dHxy[0,k]
                Bv[1,2*k  ]=dHxy[1,k]
                Bv[2,2*k+1]=dHxy[0,k]
                Bv[3,2*k+1]=dHxy[1,k]
            
            #Galerkin strain integration
            #Calculate deformation gradient Fij, for that
            #Calculate Velocity gradient Lij
            #Lij=dvi/dxj(2.4) 
            #According to B.11 d(phi)/dxj=J-1(ij) dphi/drj = Bvjk Vk
            dVxy=Bv*Ve #(4x1)(vx,x vx,y vy,x vy,y)T 
            
            #Stabilization factor tau 2.26
            #tau=beta*he/(2|v|)
            #See beta estability paramter
            #LM (2.33 p23)
            LM[0,0]=LM[2,2]=dVxy[0]
            LM[0,l]=LM[2,3]=dVxy[1]
            LM[1,0]=LM[3,2]=dVxy[2]
            LM[1,1]=LM[3,3]=dVxy[3]
            
            #BL interpolators BLijk (4,4,8) (B.17 p165)
            
            #Set tractions tP (2.34)
            
            
            #Before Strain tensors, DEFORMATIONS
            #With evp = f(sigma,s) (Sec. 2.4.2)
            #eps_vp=A sinh (psi (sigma/s))^(1/m) (2.57 p27)
            #Calculate Rate of Def Tensor D (2.13, 2.14)
            #D(th)ij=alpha vk dT/dxk deltaij 2.13
            #IN THIS EXAMPLE THIS DTH 0
    #D(vp)ij=sqrt(3/2) e. vp Nij  2.14
    #Nij Direction of plastic flow
    #
    #Calculate Leij = Lij - D(th)ij - D(vp) ij (2.10-2.12)
    #
            #Calculate La (ij)
            #Laij=F(-1)ki F(-1)kl Le lj
            #Calculate Almansi deformation gradient E (A.5)
            #Ea ij= 1/2(Lki F(-1)lk F(-1)LJ +F(-1)ki F(-1)KL Llj )

            w=0.25 #TO MODIFY
            #Calculate sigma
            #Calculate Piola Kirchoff Pi (2.31) Gij Cjk˙¯-Gij LM (jk) sig(k) + 
            #Attention double contraction
            P=G*c*E-G*LM*sig
            #print (B)
            #K+=(B.transpose()*c*B*w)
            #print (K)
    #print (K)

#Boundary conditions
#Numerated as in Bathe
for i in range(8):
    K[4,i] = K[i,4] = 0.0
    K[5,i] = K[i,5] = 0.0
    K[7,i] = K[i,7] = 0.0
	
K[4,4] = K[5,5] = K[7,7] = 1.;
R=[0,0,1000.0,0,0,0,0,0]

# #Numerated as in deal.ii
# for i in range(8):
    # K[0,i] = K[i,0] = 0.0
    # K[1,i] = K[i,1] = 0.0
    # K[3,i] = K[i,3] = 0.0
	
# K[0,0] = K[1,1] = K[3,3] = 1.;
# R=[0,0,0,0,1000.0,0,0,0]

print (K)

U=linalg.solve(K, R)
print ("Results")
print(U)


