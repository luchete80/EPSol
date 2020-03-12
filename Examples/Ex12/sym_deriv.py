from sympy import *

s,h0,A,s_,s_ast,m,n,psi,Uf,Ufvp,Dvp = symbols('s h0 A s_ s_ast m n psi Uf Ufvp Dvp')

#From evp=f(s,sig)

z = symbols('z', cls=Function)
Ee = symbols('Ee', cls=Function)
sig=z(Uf,Ufvp,Ee)
f= A*Pow(psi*sinh(s/sig),1/m)

Dvp=sqrt(3.2)*f
#1-s/s* > 0
g1 = h0*(1-s/(s_*Pow(f/A,n)))*f

print("f's: ",diff(f,s))
print("f's: ",diff(f,s))
print("-------------------------")
print("g1's: ",diff(g1,s))
print("-------------------------")
print("g1'Uf: ",diff(g1,Uf))
print("-------------------------")

from sympy import *
r, t = symbols('r t') # r (radius), t (angle theta)
f = symbols('f', cls=Function)
x = r* cos(t)
y = r* sin(t)
g = f(x,y)
#print(Derivative(g,r, 2).doit())

