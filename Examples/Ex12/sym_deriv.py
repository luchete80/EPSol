from sympy import *

s,h0,A,s_,s_ast,m,n,psi = symbols('s h0 A s_ s_ast m n psi')

#1-s/s* > 0
#From evp=f(s,sig)
sig=s*s
f= A*Pow(psi*sinh(s/sig),1/m)
g1 = h0*(1-s/(s_*Pow(f/A,n)))*f

print("f's: ",diff(f,s))
print("g1's: ",diff(g1,s))

from sympy import *
r, t = symbols('r t') # r (radius), t (angle theta)
f = symbols('f', cls=Function)
x = r* cos(t)
y = r* sin(t)
g = f(x,y)
#print(Derivative(g,r, 2).doit())

