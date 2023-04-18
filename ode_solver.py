# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a ode_solver script file.
"""


import numpy as np
import matplotlib.pyplot as plt
import math
import time


def euler(func,u0,t0,deltat):
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


#%%

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
    soln : list
        the first element in the list is the numpy array containing all solutions for each timestep
        the second element is the numpy array containing all time from t1 to t2, one stepsize interval

    """
    soln = np.reshape(x1,(-1,1))
    deltat = t2-t1
    # change the deltat until it is smaller than deltat_max
    while deltat >= deltat_max:
        deltat = deltat/10 
    num_steps = (t2-t1)/deltat
    t = np.linspace(t1,t2,int(num_steps)+1)
    
    # method defines the method to solve ode
    if method == 'euler':
        for i in range(0,int(num_steps)):
            #use euler method to solve ode (one step)
            x1,t1 = euler(func,x1,t1,deltat)
            #store solution 
            xr= np.reshape(x1,(-1,1))
            soln = np.column_stack((soln,xr))
    elif method == 'RK4':
        for i in range(0,int(num_steps)):
            x1,t1 = RK4(func,x1,t1,deltat)
            xr = np.reshape(x1,(-1,1))
            soln = np.column_stack((soln,xr))
    elif method == 'midpoint':
        for i in range(0,int(num_steps)):
            x1,t1 = midpoint(func,x1,t1,deltat)
            xr = np.reshape(x1,(-1,1))
            soln = np.column_stack((soln,xr))
    else:
        x1 = 'nan'
        method = 'nan'
        print('Wrong input method,please try again')
    
    #soln.append(str(x1))
    #soln.append(str(t2))
    #soln.append(method)
    return [x1,soln,t]


#%%


if __name__=='__main__':
    
    #define ode function
    def func(x,t):
        dxdt = x 
        return dxdt
    
    def func1(u,t):
        x,y = u
        dxdt = y
        dydt = -x
        return np.array([dxdt,dydt])

    z = solve_to(func,1,0,1,'RK4',0.001)
    print('The solution of ode at t = '+ str(z[2][-1]) + ' is ' + str(z[0]))

#%%

    #Error of RK4 and Euler method against stepsize deltat
    
    error_euler = []
    error_RK4 = []
    Time = []
    Time2 =[]

    #step size
    deltat = [1,0.5,0.2,0.1,0.05,0.02,0.01,0.005,0.002,0.001]

    for ele in deltat:
        t0,x0,y0 =0,1,1
        num_steps = 1/ele
    
        for i in range(0,int(num_steps)):
            x0,t0 = euler(func,x0,t0,ele)
            y0,t0 =RK4(func,y0,t0,ele)

        x_true = math.exp(1)
        error_euler.append(abs(x_true - x0))
        error_RK4.append(abs(x_true-y0))

    #plot the errors against deltat of two methods
    plt.plot(deltat,error_euler,deltat,error_RK4)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('stepsize')
    plt.ylabel('Error')
    plt.title('Error of RK4 and Euler method against stepsize deltat')
    plt.legend(["Euler","RK4"], loc ="lower right")





    #%%
    #output the errors,see the value of deltats give the same error

    #print(error_euler)
    #print(error_RK4)
    #%%
    #print the error,  we see that the error is the same around 0.001 when deltat is 0.01 for Euler method and 0.5 for RK4. Calculate the run time


    #calculate the run time for euler 
    x0 = (1,1)
    t0 = 0
    t1 = time.time()
    for i in range(0,int(1/0.01)):
        x0,t0 = euler(func1,x0,t0,0.01)
    Time1 = time.time()-t1
    print('The running time of Euler method solving ode with an error of 0.003 is ' + str(Time1))


    #calculate the run time for RK4
    x0 = (1,1)
    t0 = 0
    t2 = time.time()
    for i in range(0,int(1/0.5)):
        x0,t0 = RK4(func1,x0,t0,0.5)
    Time2 = time.time()-t2
    print('The running time of RK4 method solving ode with an error of 0.003 is ' + str(Time2))
        
    #%%
    
    
        
    # solution to second order ode using different stepsize
    z1 = solve_to(func1,(1,1),0,1,'euler',0.1)
    z2 = solve_to(func1,(1,1),0,1,'RK4',0.1)
    z3 = solve_to(func1,(1,1),0,1,'euler',0.01)
    z4 = solve_to(func1,(1,1),0,1,'RK4',0.01)
    
    print(z1)
    
    # plot the ode solution x against t
    plt.figure()
    plt.plot(z1[2], z1[1][0], 'b',label='Euler(deltat=0.1)')
    plt.plot(z2[2], z2[1][0], 'k',label='RK4(deltat=0.1)')
    plt.plot(z3[2], z3[1][0], 'g-',label='Euler(deltat=0.01)')
    plt.plot(z4[2], z4[1][0], 'y-',label='RK4(deltat=0.01)')
    plt.plot(z4[2],np.sin(z4[2])+np.cos(z4[2]),'r',label='real solution')
    plt.ylabel('x')
    plt.xlabel('t')
    plt.legend()
    plt.title('solution to second order ode using different stepsizes')
    plt.show()


    # plot x against t and y against t in time interval [0,10]
    z5 = solve_to(func1,(1,1),0,10,'RK4',0.01)
    plt.plot(z5[2], z5[1][0], z5[2], z5[1][1],label='RK4(deltat=0.01)')

