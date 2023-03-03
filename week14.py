# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#use euler's method to solve a ode probelm
import numpy as np
import matplotlib.pyplot as plt
import math



def euler_step(x0,y0,deltat):
    y1 = y0+deltat*x0
    return y1


#the solution to the ode equation

#initial values

error = []

#step size
deltat = [1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001]

for t in deltat:
    x0 =0
    y0 =1
    num_steps = 1/t
    i =1
    while i <= num_steps:
        y0 = euler_step(x0,y0,t)
        x0 = x0+t
        i += 1
    
    y_true = math.exp(1)
    error.append(abs(y_true - y0))

print(y0)
print(error)
plt.plot(deltat,error)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('stepsize')
plt.ylabel('Error')
