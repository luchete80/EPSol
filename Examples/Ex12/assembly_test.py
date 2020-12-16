from numpy import *

nex=2
ney=1

numel=nex*ney
vnrow=zeros(20) #Dof connectivity
vncol=zeros(20) #Dof connectivity
numvars=2 #1: Only F tensor, 2: F and internal variable s
numnodes=(nex+1)*(ney+1)
node=zeros((numnodes, 2))
vnxy=zeros((numnodes, 2))
elnodes=zeros((numel, 4))#Connectivity
elnodes.astype(int)

var_edof=zeros(2)
var_dim =[4,1]

vrowinc=0
#Assembly Matrix
for vrow in range(numvars): #Variables
    ir=0
    imax=int(var_dim[vrow])
    for n in range (4): #Nodes
        for i in range(imax): 
            d=elnodes.astype(int)[e][n]
            print("vrowinc,d,a+b",vrowinc,d,vrowinc+var_dim[vrow]*d+i)
            vnrow[ir]=vrowinc+var_dim[vrow]*d+i
            ir=ir+1
    
    print ("vrow vnrow",vrow,vnrow)
    vcolinc=0        
    for vcol in range(numvars): #Variables
        print("vcol",vcol)
        jmax=int(var_dim[vcol])
        print("imax, jmax",imax,jmax)
        #Store vn vectors
        ic=0
        for n in range (4): #Nodes 
            for j in range(jmax):
                d=elnodes.astype(int)[e][n]
                print("vcolinc",vcolinc)
                vncol[ic]=vcolinc+var_dim[vcol]*d+j
                ic=ic+1
                    
        for row in range(4*imax):
            for col in range(4*jmax):
                Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]]=  Kglob[vnrow.astype(int)[row],vncol.astype(int)[col]]+(
                                                                      Kt[vrow][vcol][row,col])
        vcolinc+=numnodes*var_dim[vcol]
    
    Rglob[vnrow.astype(int)[row]]=R[vrow][row]
    
    vrowinc+=numnodes*var_dim[vrow]