from sympy import Symbol
x, y = Symbol('x y')

f = x + y
vx= 0.1/(sqrt(x*x+y*y))*cos(atan(y/x))
dvxdy=diff(vx,y)
print("fsym",f)
print("dvxdy",dvxdy)
f.subs({x:10, y: 20})
print("f",f)
    
