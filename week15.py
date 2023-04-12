# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 18:39:32 2023

@author: YHU
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.optimize import fsolve

#define the predator_prey function as a system of odes
def predator_prey(x0,y0,a,b,d,deltat):
    dxdt = x0*(1-x0)-(a*x0*y0)/(d+x0)
    dydt = b*y0*(1-y0/x0)
    x1 = x0 + dxdt*deltat
    y1 = y0 + dydt*deltat
    return x1,y1


#solve the x and y corresponding to number of preys and predators
#plot the gragh of their population against time with different value of b
def predator_prey_plot(b):
        
    a = 1
    d = 0.1 
    deltat = 0.01
    x0,y0 = 0.05,0.05
    soln_x = [x0]
    soln_y = [y0]
    t = np.linspace(0,200,int(200/deltat)+1)
    for i in range(int(200/deltat)):
        x0,y0 = predator_prey(x0,y0,a,b,d,deltat)
        soln_x.append(x0)
        soln_y.append(y0)

    plt.figure()
    plt.plot(t,soln_x)
    plt.plot(t,soln_y)
    plt.legend(['Prey','Predator'])
    plt.yticks(np.arange(0,1,0.1))
    plt.xlabel('Time')
    plt.ylabel('Population')
    plt.show()
    return soln_x, soln_y

#solutions by using different b
soln4_x,soln4_y = predator_prey_plot(0.4)
soln3_x,soln3_y = predator_prey_plot(0.3)
soln2_x,soln2_y = predator_prey_plot(0.2)
soln1_x,soln1_y = predator_prey_plot(0.1)

# use different values of b, by looking at the plots, we see that for long-time limit
# when b <0.26, (b= 0.1/b=0.2), the x and y are oscillating and have limit cylces
# when b> 0.26 (b =0.3/b=0.4), the x and y converge together to a equilibrium point 
plt.plot(soln1_x,soln1_y,label = 'b = 0.1')
plt.plot(soln2_x,soln2_y,label = 'b = 0.2')
plt.plot(soln3_x,soln3_y,label = 'b = 0.3')
plt.plot(soln4_x,soln4_y,label = 'b = 0.4')
plt.xlabel('x')
plt.ylabel('y')
plt.legend()




# define the predator_prey ode 
def predator_prey2(x0,y0,a,b,d,deltat):
    dxdt = x0*(1-x0)-(a*x0*y0)/(d+x0)
    dydt = b*y0*(1-y0/x0)
    return [dxdt,dydt]

#x nullcline : dx/dt = 0 reference(bb)
xval = np.linspace(0.01,1,201)
yval = np.zeros(np.size(xval))
for (i,V) in enumerate(xval):
    result = fsolve(lambda Y:predator_prey2(V,Y,1,0.1,0.1,0.01)[0],1)
    #if result.success:
    yval[i] = result
    #else:
        #yval[i] = nan
plt.figure()
plt.plot(xval,yval)
print(yval)


#y nullcline  dy/dt = 0
xval = np.linspace(0.01,1,201)
yval = np.zeros(np.size(xval))
for (i,V) in enumerate(xval):
    result = fsolve(lambda Y:predator_prey2(V,Y,1,0.1,0.1,0.01)[1],1)
    #if result.success:
    yval[i] = result
    #else:
        #yval[i] = nan
plt.plot(xval,yval)
print(yval)


#find the equilibrium point where dy/dt = dx/dt = 0
x,y = fsolve(lambda u :predator_prey2(u[0],u[1],1,0.1,0.1,0.01),(0.3,0.3))
print('The equilibrium point is at '+ '('+str(x)+','+ str(y)+')')



#numerical shooting method with fsolve function. this function is specified to solve predator ode.
# to solve higher dimension system odes, we will improve it later. 
#x0,y0 is the initial guess of solution that will give limit cycles. T is the guessed period
# To have a better guess, we look at x,y againt t plots, the occilation of limit cycles starts at
#around x = 0.3 and y =0.3 with a period around 20
 
def numerical_shooting(x0,y0,T):
    a = 1
    b = 0.1
    d = 0.1
    deltat = 0.01
    p = x0
    q = y0
    dxdt = x0*(1-x0)-(a*x0*y0)/(d+x0)
    for i in range(int(T/deltat)):
        p,q=predator_prey(p,q,a,b,d,deltat)
    return [x0-p,y0-q,dxdt]
    
def solve(x,y,T):
    
    x,y,z = fsolve(lambda u:numerical_shooting(u[0],u[1],u[2]),(x,y,T))
    
    print('The initial value gives a limit cycle is x = ' + str(x) +','+ 'y = '+ str(y) +' with a period '+ str(z) )
    return np.array([x,y,z])

solve(0.3,0.3,10)

# use fsolve in scripy, the pros are this is an existing function that is well-built to find roots,
# the cons are we need give initial guess and it may not converge and it may give reuslt that is wrong.like T less than 0   
