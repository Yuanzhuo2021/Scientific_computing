# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#use euler's method to solve a ode probelm
import numpy as np
import matplotlib.pyplot as plt
import math



def euler_step(x0,t0,deltat):
    x1 = x0+deltat*t0
    return x1


#the solution to the ode equation

#initial values

error = []

#step size
deltat = [1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001]

for ele in deltat:
    t0 =0
    x0 =1
    num_steps = 1/ele
    i =1
    while i <= num_steps:
        x0 = euler_step(x0,t0,ele)
        t0 = t0+ele
        i += 1
    
    x_true = math.exp(1)
    error.append(abs(x_true - x0))

print(x0)
print(error)
plt.plot(deltat,error)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('stepsize')
plt.ylabel('Error')
plt.title('Error of Eulers method against stepsize deltat')



def solve_to(x1,t1,t2,deltat_max):
    deltat = t2-t1
    while deltat >= deltat_max:
        deltat = deltat/10
     
    num_steps = (t2-t1)/deltat
    i =1
    while i <= num_steps:
        x1 = euler_step(x1,t1,deltat)
        t1 = t1+deltat
        i += 1
    return x1
   

x2 = solve_to(1,0,1,0.01)
print(x2)
