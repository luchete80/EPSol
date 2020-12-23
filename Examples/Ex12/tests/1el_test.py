# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
#import math
from numpy import *
import numpy.matlib #for zeros

import array as arr

#******************************************************
#----------------------------
#Input Data------------------
lx=1.
ly=1
nex=1
ney=1
numit=2 

#-------------
numvars=1 #1: Only F tensor, 2: F and internal variable s

##********************** INPUT END 

#***** MESH**********
dx=lx/nex
dy=ly/ney
numel=nex*ney

vnrow=zeros(20) #Dof connectivity
vncol=zeros(20) #Dof connectivity

numnodes=(nex+1)*(ney+1)
node=zeros((numnodes, 2))
vnxy=zeros((numnodes, 2))
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

#BEFORE CONVERT TO CARTESIAN!
#Only for this example   
print("************VELOCITIES***************")
vnxy[0,0]=0.;vnxy[0,1]=1.
vnxy[1,0]=0.1;vnxy[1,1]=1.
vnxy[2,0]=0.;vnxy[2,1]=1.
vnxy[3,0]=0.1;vnxy[3,1]=1.

for n in range (numnodes):
    print("vxy ",n,":",vnxy[n,0],vnxy[n,1])
 
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
#Element dimension and DOF PER VARIABLE! 
var_edof=zeros(2)
var_dim =[4,1]
ndof=0
#Formulation type and DOFs
for i in range(numvars): 
    ndof += var_dim[i]
print ("dofs: ", ndof)
for i in range(numvars):
    var_edof[i]=4*var_dim[i]  
    
#print ("var_edof",var_edof)

edof=ndof*4
dof=ndof*numnodes
print ("GLOBAL DOF: ",dof)
#X is combines per row to calculate J
p=1.0/1.732050807568877
gauss=[-p,p]
#Numerated as in Bathe
X2=numpy.matlib.zeros((4, 2))

#Main variables values
#Nodal values
U   =matrix(numpy.matlib.zeros((44, 1)))
UV  =matrix(numpy.matlib.zeros((8, 1)))
UF  =matrix(numpy.matlib.zeros((16, 1)))
Us  =matrix(numpy.matlib.zeros((4, 1)))

#Gauss variables
#d vectors are [xx yx xy yy zz]
F       = matrix(numpy.matlib.zeros((4,1)))#Tensor form

S=matrix(numpy.matlib.zeros((4, 1)))    #Nodal internal variable
v=matrix(numpy.matlib.zeros((2, 1)))    #Gauss Point velocity

#Shape Functions
Ns=matrix(numpy.matlib.zeros((1, 4)))
Nv=matrix(numpy.matlib.zeros((2, 8)))
NsigF=matrix(numpy.matlib.zeros((4, 16)))

Dvp=matrix(numpy.matlib.zeros((4, 1)))      #From plastic deformation direction
DM =matrix(numpy.matlib.zeros((5, 5)))   #4.24
temp2x2=matrix(numpy.matlib.zeros((2, 2)))

#Derivatives
dHxy=matrix(numpy.matlib.zeros((2, 4)))
Bs=matrix(numpy.matlib.zeros((2, 8)))
Bv=matrix(numpy.matlib.zeros((4, 8)))
#BsigF=[matrix(numpy.matlib.zeros((4, 16))),matrix(numpy.matlib.zeros((4, 8)))]
BsigF=arange(128).reshape(4,16,2) #

temp4x16=matrix(numpy.matlib.zeros((4, 16)))

B4i=arange(32).reshape(4,4,2) #
#(4,16,2)
#print(BsigF[0])
LM =matrix(numpy.matlib.zeros((4, 4)))


R   =[  matrix(numpy.matlib.zeros((16, 1))),
        matrix(numpy.matlib.zeros(( 4, 1)))]
RF  =matrix(numpy.matlib.zeros((16, 1)))
Rsig=matrix(numpy.matlib.zeros((16, 1)))
Rs  =matrix(numpy.matlib.zeros((4, 1)))
Rv  =matrix(numpy.matlib.zeros((8, 1)))


#These are the same but reorganized
dVxy=zeros(4)
BL      = arange(128).reshape(4,4,8)            #Eqns 2.33, B.17
temp8x1 = matrix(numpy.matlib.zeros((8, 1))) 


P=matrix(numpy.matlib.zeros((4, 1))) 

#Global matrices
#Uglob=matrix(numpy.matlib.zeros((ndof*numnodes, ndof*numnodes)))
Kglob=matrix(numpy.matlib.zeros((dof, dof))) 
dUglob=zeros(dof)
Uglob=zeros(dof)
Rglob=zeros(dof)

class bMatrix: #Block matrix
    
    def __init__(self,i,j):
        m=matrix(numpy.matlib.zeros((i, j)))
 
Kt=[
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[0], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[0], var_edof.astype(int)[1])))]
     ,
     [matrix(numpy.matlib.zeros(( var_edof.astype(int)[1], var_edof.astype(int)[0]))), matrix(numpy.matlib.zeros((var_edof.astype(int)[1], var_edof.astype(int)[1])))]
    ] 
#---------------------------------------------------
#Before start, assign Uglob values
#---------------------------------------------------
#To MODIFY put in a function     

for n in range(numnodes):
    #Velocities 
    iF=ndof*n
    #Initial deformation gradients as identity??
    Uglob[iF  ]=Uglob[iF+3]=1
    Uglob[iF+1]=Uglob[iF+2]=0#xy and yx        

                        
print ("Initial Uglob", Uglob)

#-------------------------------------------
it=0

## ------------------------------------------
## Newton Rhapson Loop
## ------------------------------------------
while (it < numit):
    
#ELEMENT LOOP  ----------------
    #for e in range (1):
    for e in range (numel): # TO MODIFY 
        #Obtain Ve from global
        Kel=0.
        for n in range(4):
            X2[n]=node[elnodes.astype(int)[e][n]]
        #print ("Element ", e)
        #print ("Element Nodes")
        #print (X2)
        
        for ig in range(2):
            for jg in range(2):
                rg=gauss[ig]
                sg=gauss[jg]

                #Numerated as in Bathe
                Ns  =0.25*matrix([(1+sg)*(1+rg),(1-rg)*(1+sg),(1-sg)*(1-rg),(1-sg)*(1+rg)])   
                dHrs=matrix([[(1+sg),-(1+sg),-(1-sg),(1-sg)], [(1+rg),(1-rg),-(1-rg),-(1+rg)] ])
                #print ("Ns",Ns)
                #Numerated as in deal.ii
                #dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
                dHrs/=4.
                J=dHrs*X2
                dHxy=linalg.inv(J)*dHrs
                detJ=linalg.det(J)
                #Calculate shape functions
                #Bs=J-1 dHrs(B.13)
                Bs=dHxy
                for k in range(4):
                    #shape functions
                    Nv[0,2*k  ]=Nv[1,2*k+1]=Ns[0,k]
                    for j in range(4):
                        NsigF[j,4*k+j]=Ns[0,k]

                    #derivatives Bv (B.14)
                    Bv[0,2*k  ]=dHxy[0,k]
                    Bv[1,2*k  ]=dHxy[1,k]
                    Bv[2,2*k+1]=dHxy[0,k]
                    Bv[3,2*k+1]=dHxy[1,k]
                
                #print("NsigF",NsigF)
                ################            
                ## 
                for i in range(4):
                    for l in range(4):
                        for m in range(4):  
                            for n in range(2):
                                if (l==m):
                                    B4i[l,m,n]=Bs[n,i]
                                else:
                                    B4i[l,m,n]=0.              
                    for l in range(4):
                        for m in range(4):  
                            for n in range(2):
                                BsigF[l,4*i+m,n]=B4i[l,m,n]
      
                                
                #Interpolate velocity
                #INCREMENT GLOBAL VELOCITY FROM INCREMENTS!!!
                #CHANGE F TO ASSEMBLE IN THE SAME PLACE FOR BOTH FORMS
                juf=0
                iv=0
                for n in range (4): #Element nodes
                    gn=elnodes.astype(int)[e][n]#globnod
                    
                    for i in range (2):    #Velocity is var 0
                        #print("UV,i+iv,node",i+iv,gn)
                        UV[i+iv,0]=vnxy[gn,i]
                        
                    for j in range (var_dim[0]):
                        UF  [j+juf,0]=Uglob[ndof*gn+j]
                        #print("UF(j+juf,dof)",j+juf,ndof*gn+j)
                    juf+=var_dim[0]
                    iv=iv+2

                print("UV",UV)
                print("UF",UF)
                
                v  =Nv*UV #[2x8 x (8x1)]
                print ("v",v)
                s  =float(Ns*Us)
                F  =NsigF*UF #[(4x16)*(16x1) =(4x1)]
                
                print ("F",F)
                 
                #Galerkin strain integration
                #Calculate deformation gradient Fij, for that
                #Calculate Velocity gradient Lij
                #Lij=dvi/dxj(2.4) 
                #According to B.11 d(phi)/dxj=J-1(ij) dphi/drj = Bvjk Vk
                
                #Formulation 1
                #Calculate La (ij) 2.9 
                #Laij=F(-1)ki F(-1)kl Le lj
                
                #Calculate Leij = Lij - D(th)ij - D(vp) ij (2.10-2.12)
                #Calculate Almansi deformation gradient E (A.5)
                #Calculate 
                #Ea ij= 1/2(Lki F(-1)lk F(-1)LJ +F(-1)ki F(-1)KL Llj )
                dVxy=Bv*UV #(4x8)*(8x1)=(4x1) (vx,x vx,y vy,x vy,y)T 
                
                print ("dVxy",dVxy)
                
                #Stabilization factor tau 2.26
                beta=1.
                he=(lx+ly)/2. #ONLY FOR THIS EXAMPLE
                tau=float(beta*he/(2.*sqrt(v[0]*v[0]+v[1]*[1])))
                #tau=1.
                print("tau",tau)
                #See beta estability paramter
                #LM (2.33 & 4.23 p23 & p91)
                #Attention: SE DIFFERENCES WITH L_ in 2.28
                LM[0,0]=LM[2,2]=dVxy[0]
                LM[0,1]=LM[2,3]=dVxy[1]
                LM[1,0]=LM[3,2]=dVxy[2]
                LM[1,1]=LM[3,3]=dVxy[3]
                
                    
                w=1. #TO MODIFY
                #Calculate sigma
                #2.31 Comes first from 2.2 (Then 2.31)
                #Pij=vk d(sig ij)/xk - dvi/dxk sig(kj) + dvk/dxk sig ij
                #Calculate Piola Kirchoff Pi (2.31) Gij Cjk-Gij LM (jk) sig(k)+
                #Attention double contraction
                #P=G*c*E-G*LM*sig+(LM[0,0]+LM[1,1])*G*sig

                wJ=w*detJ
                
             
                visc=1.
       
                # *****************************************************************
                #RESIDUALS ******************* 2.26 to 2.39 *****************************
                #Construct vk Bsig mik
                for m in range(4):
                    for i in range(16):
                            temp4x16[m,i]=0.
                for m in range(4):
                    for i in range(16):
                        for k in range(2):
                            temp4x16[m,i]=temp4x16[m,i]+BsigF[m,i,k]*v[k,0]
               
                
                R[0]   =(NsigF+temp4x16*tau).transpose()*(temp4x16*UF-LM*NsigF*UF)*wJ
              
                #R Assembly            
                #TANGENT MATRIX   
                #PAGES 25 y 94
                
                #----------------------------------------------------                   
                               
                #dRFdUF  4.36

                Kt[0][0]= Kt[0][0]+(
                                (NsigF+temp4x16*tau).transpose()*
                                (temp4x16-LM*NsigF)
                                )*wJ
                               
        
        #print ("Nv",Nv)
        print("Kt[0][0]",Kt[0][0]) 
        
        vrowinc=0
        #Assembly Matrix
        for vrow in range(numvars): #Variables
            ir=0
            imax=int(var_dim[vrow])
            for n in range (4): #Nodes
                for i in range(imax): 
                    d=elnodes.astype(int)[e][n]
                    #print("vrowinc,d,a+b",vrowinc,d,vrowinc+var_dim[vrow]*d+i)
                    vnrow[ir]=vrowinc+var_dim[vrow]*d+i
                    ir=ir+1
            
            #print ("vrow vnrow",vrow,vnrow)
            vcolinc=0        
            for vcol in range(numvars): #Variables
                
                jmax=int(var_dim[vcol])
                #print("imax, jmax",imax,jmax)
                #Store vn vectors
                ic=0
                for n in range (4): #Nodes 
                    for j in range(jmax):
                        d=elnodes.astype(int)[e][n]
                        #print("vcolinc",vcolinc)
                        vncol[ic]=vcolinc+var_dim[vcol]*d+j
                        ic=ic+1
                
                #print("vcol vncol",vcol,vncol)
                for row in range(4*imax):
                    Rglob[vnrow.astype(int)[row]]+=R[vrow][row]
                    for col in range(4*jmax):
                        #print("vnrow(row)vncol(col)",vnrow[row],vncol[col]) 
                        Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]] =  Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]]+(
                                                                              Kt[vrow][vcol][row,col])
                vcolinc+=numnodes*var_dim[vcol]
            
            
            
            vrowinc+=numnodes*var_dim[vrow]
        

    ##---------------------------------------------------------------
    ##Boundary conditions
    ##---------------------------------------------------------------
    #Velocity DOFs 
    #In this example velocities are known
    #AT INLET(left nodes):
    # F=I , sigma = 0
    dnode=(nex+1)    
    for dy in range(ney+1): 
        inode=dy*dnode
        #print("node",inode)   
        #Deformation gradient F
        for i in range ( var_dim [ 0 ] ):
            idof = var_dim[0] * inode + i
            print ("idof",idof)
            for j in range(dof):
                Kglob[ idof , j ] = 0
                #Rglob[ j ] -= Kglob[j,idof] * 0 #dU=0, U=1(idof)
                Kglob[ j ,idof ] = 0
            

            Kglob[idof,idof] = 1
            Rglob[idof  ] = 0           #F INCREMENT (dF) IS NULL!!!!!
                      

    print("KGLOB\n")
    for i in range (dof):
        for j in range (dof):
            print(Kglob[i,j], end = " ")
            
        print("\n")
    print("Rglob",Rglob)


#print (K)
    for i in range (dof):
        Uglob[i]=Uglob[i]+dUglob[i]
    
    dUglob=linalg.solve(Kglob, Rglob)
    #print("it %d, dUglob",it, dUglob)  
    #print("Uglob", Uglob)
    
    
    #TOTAL BOUNDARY CONDITIONS FOR UF and s calculations
    dnode=(nex+1)    
    for dy in range(ney+1): 
        inode=dy*dnode
        #Deformation gradient F
        idof = var_dim[0] * inode 
        print ("inode, idof",inode, idof)
        Uglob[ idof     ] = Uglob[ idof + 3 ] = 1
        Uglob[ idof + 1 ] = Uglob[ idof + 2 ] = 0

        
        #Sigma is zero, Internal variable s ,      
        if numvars == 2:
            idofs = idof + var_dim[0]
            Uglob[ idofs ] = mat_s0
 
    
    it+=1
    print ("Iteration: ",it, "from ", numit)
    


    
print ("Results")
print("Uglob", Uglob)

file= open("output.vtu","w+")
file.write("<?xml version=\"1.0\"?>\n")
file.write("<VTKFile type=\"UnstructuredGrid\">\n")
file.write("<UnstructuredGrid>\n")
file.write("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n" %(numnodes,numel) )
file.write("<Points>\n")
file.write("<DataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\" >\n")
for i in range (numnodes):
    file.write("%f\n %f 0.\n" %(node[i,0],node[i,1]))
file.write("</DataArray>\n")
file.write("</Points>\n")
file.write("<Cells>\n")
file.write("<DataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\" >\n")
for e in range (numel):
    for n in range (4):
        file.write("%d " %(elnodes.astype(int)[e][n]))
    file.write("\n")
file.write("</DataArray>\n")
file.write("<DataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\" >\n")
offset=4
for e in range (numel):
    file.write("%d " %(offset))
    offset+=4
file.write("\n</DataArray>\n")
file.write("<DataArray Name=\"types\" type=\"UInt8\" format=\"ascii\" >\n")
for e in range (numel):
    file.write("10 ")
file.write("\n</DataArray>\n");
file.write("</Cells>\n")
file.write("<PointData Scalars=\"scalars\" format=\"ascii\">\n")
file.write("<DataArray Name=\"DefGrad_F\" NumberOfComponents=\"%d\" type=\"Float32\" format=\"ascii\" >\n" %(ndof))
v=0
for n in range (numnodes):
    for d in range (ndof):
        print("v,Uglob[v]",v,Uglob[v])
        file.write("%f " %(Uglob[v]))
        v=v+1
    file.write("\n")
file.write("\n</DataArray>\n")
file.write("<DataArray Name=\"Vel\" NumberOfComponents=\"2\" type=\"Float32\" format=\"ascii\" >\n")
v=0
for n in range (numnodes):
    for d in range (2):
        file.write("%f " %(vnxy[n,d]))
    file.write("\n")
    #i+=var_dim[0]
file.write("\n</DataArray>\n")
file.write("</PointData>\n")
file.write("</Piece>\n")
file.write("</UnstructuredGrid>\n")
file.write("</VTKFile>\n")
file.close
