# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 04:51:56 2023

@author: YHU
"""

import numpy as np
from scipy.optimize import fsolve,root
import matplotlib.pyplot as plt


def algebraic_cubic(u):
    x,c = u
    return x**3 - x + c

def continuation(func,rang,guess):
    S = 1000
    step = (rang[1]-rang[0])/S
    C = []
    solution = []
    for i in range(S):
        c = rang[0] + step*(i)
        result = root(lambda x: func((x,c)),guess)
        if result.success:
            C.append(c)
            solution = np.hstack((solution,result.x))
            result = guess
        else:
            print(result.message)
            print('Converge fails at c = ',c)
            break
    return [solution,C]

z = continuation(algebraic_cubic,(-2,2),2)

plt.plot(z[1],z[0])
plt.xlabel('Parameter')
plt.ylabel('Solution')


#%%



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
       

#%%


# delta is the secant of state vector
#delta = np.array([z[0][1],z[1][1]])-np.array([z[0][0],z[1][0]])



def psuedo_arclength_continuation(func,rang,guess):
    
    
    z = continuation(func,rang,guess)
    if len(z[0])>=2 and len(z[1])>=2:
        # get the first two solns using natural parameter continuation
        u0 = np.array([z[0][0],z[1][0]])
        u1 = np.array([z[0][1],z[1][1]])  # solution, C
        
        
       # print(delta) 
        soln = [z[0][0],z[0][1]]
        param = [z[1][0],z[1][1]]
        
        def pseudo_equation(u,u2):
            return [np.dot((u-u2),delta),func(u)]
        
        for i in range(1000):
            
            delta = (u1-u0)
            u2 = u1 + delta # estimated value
            true_value = root(lambda u: pseudo_equation(u,u2),u2) # true value of u2

            if true_value.success:
                soln.append(true_value.x[0])
                param.append(true_value.x[-1])
                u0 = u1
                u1 = true_value.x
                
            else:
                print('converge fails at c =',u1[-1])
                break
    else:
        print('Converges fail, please choose a new guess')
    return [param,soln]

zz = psuedo_arclength_continuation(algebraic_cubic,(-2,2),2)
print(zz[0])

plt.plot(zz[0],zz[1])
