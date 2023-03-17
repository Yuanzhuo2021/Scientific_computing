# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import time




#define ode function
def f(x,t):
    return x



# single step of euler's method 
def euler_step(x0,t0,deltat):
    x1 = x0+deltat*f(x0,t0)
    return x1


##single step of 4th Runge-Kutta method


    
def RK4(x0,t0,h):
    k1 = f(x0,t0) 
    k2 = f(x0+h*k1*0.5,t0+h/2)
    k3 = f(x0+h*k2*0.5,t0+h/2)
    k4 = f(x0+h,t0+h*k3)
        
    #update x0 value
    x1 = x0+h*(k1+2*k2+2*k3+k4)/6
    return x1


    


#the solution to the ode equation

#initial values

error_euler = []
error_RK4 = []
Time = []
Time2 =[]

#step size
deltat = [1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001]

for ele in deltat:
    t0 =0
    x0 =1
    y0 =1
    num_steps = 1/ele
    i =1
    
    while i <= num_steps:
        x0 = euler_step(x0,t0,ele)
        y0 =RK4(y0,t0,ele)
        t0 = t0+ele
        i += 1
        
    x_true = math.exp(1)
    error_euler.append(abs(x_true - x0))
    error_RK4.append(abs(x_true-y0))

print(x0)
print(y0)
plt.plot(deltat,error_euler,deltat,error_RK4)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('stepsize')
plt.ylabel('Error')
plt.title('Error of RK4 and Euler method against stepsize deltat')
plt.legend(["Euler","RK4"], loc ="lower right")



def solve_to(x1,t1,t2,deltat_max):
    print('Pleaswchoose a method to solve ode problem:Type 1 for RK4; Type 2 for Euler method')
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
   

z = solve_to(1,0,1,0.01)
print(z)



i = 1
x0 = 1
t0 = 0
t = time.time()
while i <= num_steps:
    x0 = euler_step(x0,t0,ele)
    t0 = t0+ele
    i += 1
Time = time.time()-t
    

