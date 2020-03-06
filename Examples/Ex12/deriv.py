from numpy import *
import numpy.matlib #for zeros

def calc_dEdU(F,Fet,Fvp):
    Fet_inv=linalg.inv(Fet)
    dEdU=[matrix(numpy.matlib.zeros((4, 16))),matrix(numpy.matlib.zeros((4, 20)))]
    #Eqn 1, comes from 
    Fet[0][0]=0.5-0.5*(1.)
    
    TF=matrix([[0,0,0,1],[0,-1,0,0],[0,0,-1,0],[1,0,0,0]])
    dEdFe_inv=matrix([[0,0,0,1],[0,-1,0,0],[0,0,-1,0],[1,0,0,0]])
    #Earr_e=
    #Earr_e=TLa*
    
    return dEdU