# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 21:59:11 2023

@author: YHU
"""
import numpy as np


# convert ode BVP into a system of algebraic equations
a = 0 
b = 10
gamma1 = 1
gamma2 = 2
N =1000
u_a = gamma1
u_b = gamma2
dx = (b-a)/N


# second order differentiation: (u[i+1] -2u[i] + u[i-1])/(dx)^2 +q[i]= 0 .Because the source equation is 0, q[i] = 0.

# use numpy to solve the system of equations. There are N-1 equations as there are N-1 unknowns.
# build A-DD
A = np.zeros((N-1,N-1))

# the first row and the last row are quite different.We deal with them seperately
A[0,0]=-2
A[0,1]=1
A[N-2,N-2]=-2
A[N-2,N-3]=1

# for the rest of the row, coefficient for u[i+1] is 1, u[i] is -2,u[i-1] is 1
for i in range(N-1-2):
    A[i+1,i]=1
    A[i+1,i+1]=-2
    A[i+1,i+2]=1
print(A)

#right side of the equations:b-DD
B = np.zeros((N-1,1))
#put initail and end conditions
B[0,0]=u_a
B[N-2,0]=u_b
print(B)

# output the soln
u = np.linalg.solve(A, -B)
x = np.linspace(a, b, N-1)


#explicite solution

def real_soln(x):
    soln = ((gamma2 - gamma1)/(b-a))*(x-a)+gamma1
    return soln

#plot 
plt.plot(x,real_soln(x),x,u)
plt.xlabel('x')
plt.ylabel('u')
plt.legend(['real solution','Finite difference'])