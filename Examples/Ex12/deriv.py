from numpy import *
import numpy.matlib #for zeros

#This function only computes dEdUF and dEdUFvp, since the other are NULL
def calc_dEdU(Fd,Fvpd,NsigF):
    #Fet_inv=linalg.inv(Fet)
    dEdU=[matrix(numpy.matlib.zeros((4, 16))),matrix(numpy.matlib.zeros((4, 20)))]
    #Eqn 1, comes from 
    #Fet[0][0]=0.5-0.5*(1.)
    TF=matrix([[0,0,0,1],[0,-1,0,0],[0,0,-1,0],[1,0,0,0]])
    dEdFed_inv=matrix(numpy.matlib.zeros((4, 5)))    #E.3
    
    dFMdUF =arange(256).reshape(4,4,16)

    Fet  = matrix(numpy.matlib.zeros((3,3)))#Tensor form    
    Fed=matrix(numpy.matlib.zeros((5, 1)))
    
    F4ed    =matrix(numpy.matlib.zeros((4, 1)))
    F4ed_inv=matrix(numpy.matlib.zeros((4, 1)))
    Fed_inv =matrix(numpy.matlib.zeros((5, 1)))
    
    Fvpt     =matrix(numpy.matlib.zeros((3, 3)))
    Fvpt_inv =matrix(numpy.matlib.zeros((3, 3)))
    Fvpd_inv =matrix(numpy.matlib.zeros((4, 1)))
    
    F4vpd_inv=matrix(numpy.matlib.zeros((4, 1)))
    
    dFeinv_dFe=matrix(numpy.matlib.zeros((5, 5)))
    dFedUF =matrix(numpy.matlib.zeros((5, 16)))
    dF4edUF=matrix(numpy.matlib.zeros((4, 16)))
    
    ddetFe_dF4ed=matrix(numpy.matlib.zeros((1, 4)))
    temp4=matrix(numpy.matlib.zeros((4, 1)))
    
    #Check with inv
    #F=[Fxx]
    print("Fd0",Fd[0,0])
    FM=matrix([[Fd[0,0],Fd[2,0],0    ,0   ],
               [Fd[1,0],Fd[3,0],0    ,0   ],
               [0    ,    0,Fd[0,0],Fd[2,0]],
               [0    ,    0,Fd[1,0],Fd[3,0]]])     
    
    
    
    # Ft[0,0]=Fd[0]
    # Ft[1,0]=Fd[1]
    # Ft[0,1]=Fd[2]
    # Ft[1,1]=Fd[3]
    # Ft[2,2]=1.
     
    Fvpt[0,0]=Fvpd[0]
    Fvpt[1,0]=Fvpd[1]
    Fvpt[0,1]=Fvpd[2]
    Fvpt[1,1]=Fvpd[3] 
    Fvpt[2,2]=Fvpd[4] 
    print ("Fvpt",Fvpt)
    
    Fvpt_inv=linalg.inv(Fvpt) 
    F4vpd_inv[0]=Fvpt_inv[0,0]
    F4vpd_inv[1]=Fvpt_inv[1,0]
    F4vpd_inv[2]=Fvpt_inv[0,1]
    F4vpd_inv[3]=Fvpt_inv[1,1]
    #Fet=Ft*Fvpt_inv #Remains thermal part
    #print(Fet)
    #E.10 
    #Fe in fact is Fe=F*Fvp-1
    F4ed=FM*Fvpd_inv         
    #Earr_e=
    #Earr_e=TLa*
    #E.9 
    for i in range (4):
        Fed[i]=F4ed[i]    
    Fed[4]=Fvpt_inv[2,2]
 
        
                
    det_Fed=Fed[0]*Fed[3]*-Fed[1]*Fed[2] #E.6
    
    #e.5
    F4ed_inv=float(-1./det_Fed)*TF*F4ed
    
    for i in range(4):
        Fed_inv[i]=F4ed_inv[i]
    
    
    #Remember E_=[Exx Eyy Ezz Exy]
    #E.3 Pag 189
    dEdFed_inv[0,0]=dEdFed_inv[3,2]=Fed_inv[0]
    dEdFed_inv[0,1]=dEdFed_inv[3,4]=Fed_inv[1]
    dEdFed_inv[1,2]=dEdFed_inv[3,0]=Fed_inv[2]
    dEdFed_inv[1,3]=dEdFed_inv[3,1]=Fed_inv[3]
    dEdFed_inv[2,4]=Fed_inv[4]

    #E.14
    for k in range(8):
        dFMdUF[0,0,k]=dFMdUF[2,2,k]=NsigF[0,k]
        dFMdUF[0,1,k]=dFMdUF[2,3,k]=NsigF[2,k] 
        dFMdUF[1,0,k]=dFMdUF[3,2,k]=NsigF[1,k]
        dFMdUF[1,1,k]=dFMdUF[3,3,k]=NsigF[3,k] 
        
    #Third Term dFedUF
    #E.12 & E.13
    for i in range(4):
        for j in range(16):
            for k in range (4):
                dF4edUF[i,j]=dFMdUF[i,k,j]*F4vpd_inv[k]
    
    #E.8
    for j in range(16):
        for i in range(4):
            dFedUF[i,j]=dF4edUF[i,j]
            dFedUF[4,j]=0
    
    #E.7
    for i,k in range(4,4):
        temp4[i,0]=temp4[i,0]+TF[i,k]*F4ed[k,0]
        
    ddetFe_dF4ed=float(-1./(det_Fed*det_Fed))*temp4+float(1./(det_Fed))*TF
    #F derivative
    #Eqn E.2
    dEdU[0]=dEdFed_inv*dFeinv_dFe*dFedUF
    
    #-----------------------------------------------------------------------------
    #Viscoplastic Fvp derivative
    #E.14
       
    return dEdU