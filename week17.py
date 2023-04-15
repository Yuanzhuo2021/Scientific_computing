# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 04:51:56 2023

@author: YHU
"""

import numpy as np
from scipy.optimize import fsolve,root
import matplotlib.pyplot as plt

def algebraic_cubic(x,c):
    return x**2 - x + c

def continuation(func,rang):
    S = 1000
    step = (rang[1]-rang[0])/S
    # first generate the first value of solution
    result= root(lambda u: func(u,rang[0]),-2)
    print(rang[0])
    print(result.x)
    solution = result.x
    C = [rang[0]]
    for i in range(S):
        c = rang[0] + step*(i+1)
        result = root(lambda u: func(u,c),result.x)
        if result.success:
            result = result
            C.append(c)
            solution = np.hstack((solution,result.x))
        else:
            print(result.message)
            print('Converge fails at c = ',c)
            break
        
    return [C,solution]

z = continuation(algebraic_cubic,(-2,2))



plt.plot(z[0],z[1])
plt.xlabel('Parameter')
plt.ylabel('Solution')






S =1000
deltat = 0.01
u1 = []
u2 = []
u0 = (1,1)
beta = 0

def Hopf_bifurcation(u,beta):
    u1,u2 = u
    du1dt = beta*u1-u2 - u1*(u1**2+u2**2)
    du2dt = u1+beta*u2 - u2*(u1**2+u2**2)
    return np.array([du1dt,du2dt])

def euler_step(u0,beta,deltat):
    u1 = u0+deltat*Hopf_bifurcation(u0,beta)
    beta = beta + 2/S
    return u1,beta

for i in range(S):
    u0,t0 = euler_step(u0,beta,deltat)
    u1.append(u0[0])
    u2.append(u0[1])

t = np.linspace(0,1,S)
beta = np.linspace(0,2,S)
print(t)
print(u1)
plt.figure()
plt.plot(beta,u1,beta,u2)
        




