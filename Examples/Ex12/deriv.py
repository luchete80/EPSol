from numpy import *
import numpy.matlib #for zeros

def calc_dEdU(Fet):
    Fet_inv=linalg.inv(Fet)
    dEdU=[matrix(numpy.matlib.zeros((4, 4))),matrix(numpy.matlib.zeros((4, 5)))]
    #Eqn 1, comes from 
    Fet[0][0]=0.5-0.5*(1.)
    
    TF=matrix([[0,0,0,1],[0,-1,0,0],[0,0,-1,0],[1,0,0,0]])
    #Earr_e=
    #Earr_e=TLa*
    
    return dEdU