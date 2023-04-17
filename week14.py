# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import time




# , the state initial values of variables in ode； t0 is the initial time； deltat is the stepsize
def euler_step(func,u0,t0,deltat):
    """
    single step of euler's method works for any dimension of ode
    
    Parameters
    ----------
    func: ode function
        An n-th order ordinary differential equation we are solving. It has to be written 
        in a system of n first-order odes.Also,the parameters should be included 
        in the 'func'.The func returns a numpy array with u1,u2,u3...un
        
        Example: A predator_prey ode function can be written as below:
        def predator_prey(u,t):
            x,y = u
            a = a
            b = b
            d = d
            dxdt = x*(1-x)-(a*x*y)/(d+x)
            dydt = b*y*(1-y/x)
            return np.array([dxdt,dydt])
        
    u0 :numpy.array 
        The state initial values of variables in ode
    t0 : int/float
        Initial time
    deltat : float
        The stepsize

    Returns
    -------
    u1 : numpy.array
        The solution to the system of odes 
    t1 : int/float
        updated t0 after one step euler 
        
    """
    u1 = u0+deltat*func(u0,t0)
    t1 = t0 + deltat
    return u1,t1


#single step of 4th Runge-Kutta method
def RK4(func,u0,t0,deltat):
    """
    single step of fourth Runge-Kutta method works for any dimension of ode
    
    Parameters
    ----------
    func: ode function
        An n-th order ordinary differential equation we are solving. It has to be written 
        in a system of n first-order odes.Also,the parameters should be included in 
        the 'func'.The func returns a numpy array with u1,u2,u3...un
        
        Example: A predator_prey ode function can be written as below:
        def predator_prey(u,t):
            x,y = u
            a = a
            b = b
            d = d
            dxdt = x*(1-x)-(a*x*y)/(d+x)
            dydt = b*y*(1-y/x)
            return np.array([du1dt,du2dt])
    
    u0 : numpy.array 
        The state initial values of variables in ode
    t0 : int/float
        Initial time
    deltat : float
        The stepsize

    Returns
    -------
    u1 : numpy.array
        The solution to the system of odes after one step of RK4
    t1 : int/float
        updated t0 after one step euler 
        
    """

    k1 = func(u0,t0)  
    k2 = func(u0+deltat*k1*0.5,t0+deltat/2)
    k3 = func(u0+deltat*k2*0.5,t0+deltat/2)
    k4 = func(u0+deltat*k3,t0+deltat)
    #update u0,t0 value
    u1 = u0+deltat*(k1+2*k2+2*k3+k4)/6
    t1 = t0 + deltat
    return u1,t1
    

# implement midpoint method
def midpoint(u0,t0,deltat):
    """
    single step of midpoint method works for any dimension of ode
    
    Parameters
    ----------
    func: ode function
        An n-th order ordinary differential equation we are solving. It has to be 
        written in a system of n first-order odes.Also,the parameters should
        be included in the 'func'.The func returns a numpy array with u1,u2,u3...un
        
        Example: A predator_prey ode function can be written as below:
        def predator_prey(u,t):
            x,y = u
            a = a
            b = b
            d = d
            dxdt = x*(1-x)-(a*x*y)/(d+x)
            dydt = b*y*(1-y/x)
            return np.array([du1dt,du2dt])
    
    u0 : numpy.array 
        The state initial values of variables in ode
    t0 : int/float
        Initial time
    deltat : float
        The stepsize

    Returns
    -------
    u1 : numpy.array
        The solution to the system of odes after one step of midpoint method
    t1 : int/float
        updated t0 after one step midpoint method

   """
    t1 = t0+0.5*deltat
    u1 = u0 + 0.5*deltat*func(u0,t0)
    k = func(u1,t1)
    t2 = t0 + deltat
    u2 = u0 + deltat*k
    return u2,t2








#define ode function
def f(x,t):
    return x

    


#the solution to the ode equation

#initial values

error_euler = []
error_RK4 = []
Time = []
Time2 =[]

#step size
deltat = [1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001,0.0005,0.0002,0.0001]

for ele in deltat:
    t0,x0,y0 =0,1,1
    num_steps = 1/ele
    i =1
    
    while i <= num_steps:
        x0 = euler_step(f,x0,t0,ele)
        y0 =RK4(f,y0,t0,ele)
        t0 = t0+ele
        i += 1
        
    x_true = math.exp(1)
    error_euler.append(abs(x_true - x0))
    error_RK4.append(abs(x_true-y0))

print(x0)
print(y0)

#plot the errors against deltat of two methods
plt.plot(deltat,error_euler,deltat,error_RK4)
plt.yscale('log')
plt.xscale('log')
plt.xlabel('stepsize')
plt.ylabel('Error')
plt.title('Error of RK4 and Euler method against stepsize deltat')
plt.legend(["Euler","RK4"], loc ="lower right")



def solve_to(func,x1,t1,t2,method,deltat_max):
    """
    This function is allow user to choose any method(euler,RK4,midpoint) to solve ode at t = t2

    Parameters
    ----------
    x1 : numpy.array
        Initial values of state
    t1 : int/float
        Initial time
    t2 : int/float
        user wants the solution at t = t2
    method : string
        user input, type 'euler' for euler's method, 'RK4' for 4th Runge Kutta method, 
        'midpoint' for midpoint method
    deltat_max : int/float
        The maximum stepsize allowed

    Returns
    -------
    soln : numpy.array
        solution to the ode at t= t2

    """
    soln = []
    deltat = t2-t1
    while deltat >= deltat_max:
        deltat = deltat/10
     
    num_steps = (t2-t1)/deltat
    if method == 'euler':
        for i in range(0,num_steps):
            x1 = euler_step(func,x1,t1,deltat)
            t1 = t1+deltat
    elif method == 'RK4':
        for i in range(0,num_steps):
            x1 = RK4(func,x1,t1,deltat)
            t1 = t1+deltat
    elif method == 'midpoint':
        for i in range(0,num_steps):
            x1 = midpoint(func,x1,t1,deltat)
            t1 = t1+deltat

    else:
            print('Wrong input,please try again')
    
                    
    soln.append(method)
    soln.append(str(x1))
    return soln

z = solve_to(1,0,1,0.001)
print('The x2 value using ' + z[0] +' method is '+ z[1])



#output the errors,see the value of deltats give the same error

print(error_euler)
print(error_RK4)

#print the error,  we see that the error is the same around 0.003 when deltat is 0.0002 for Euler method and 0.002 for RK4. Calculate the run time


i = 1
x0 = 1
t0 = 0
t1 = time.time()
while i <= 1/0.0002:
    x0 = euler_step(x0,t0,0.0002)
    t0 = t0+0.0002
    i += 1
Time1 = time.time()-t1
print('The running time of Euler method solving ode with an error of 0.003 is ' + str(Time1))


i = 1
x0 = 1
t0 = 0
t2 = time.time()
while i <= 1/0.002:
    x0 = RK4(x0,t0,0.002)
    t0 = t0+0.002
    i += 1
Time2 = time.time()-t2


print('The running time of RK4 method solving ode with an error of 0.003 is ' + str(Time2))
    


def solve_ode(deltat):
    
    x0,y0 = 1,1
    X0,Y0 = 1,1
    solution_euler = [x0]
    solution_RK4 = [X0]
    t = np.linspace(0, 1, int(1/deltat) + 1)
    
    for i in range(int(1/deltat)):
        y0,x0 = second_order_euler(y0,x0,deltat)
        Y0,X0 = second_order_RK4(Y0,X0,deltat)
        solution_euler.append(x0)
        solution_RK4.append(X0)
    return t,solution_euler,solution_RK4

t1,soln1_euler,soln1_RK4 = solve_ode(0.1)
t2,soln2_euler,soln2_RK4 = solve_ode(0.01)

plt.plot(t1, soln1_euler, 'b',label='Euler(deltat=0.1)')
plt.plot(t1, soln1_RK4, 'k',label='RK4(deltat=0.1)')
plt.plot(t2, soln2_euler, 'g-',label='Euler(deltat=0.01)')
plt.plot(t2, soln2_RK4, 'y-',label='RK4(deltat=0.01)')
plt.plot(t2,np.sin(t2)+np.cos(t2),'r',label='real solution')
plt.ylabel('x')
plt.xlabel('t')
plt.legend()
plt.title('solution to second order ode using different size of timestep')
plt.show()




deltat = 0.001
x0,t0  = 1,0
for i in range (int(1/deltat)):
    x0,t0 = midpoint(x0,t0,deltat)

print(x0)


