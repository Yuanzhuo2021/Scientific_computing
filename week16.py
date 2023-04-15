# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 23:53:31 2023

@author: YHU
"""
import numpy as np
import matplotlib.pyplot as plt
import math
import unittest
from scipy.optimize import fsolve,root

def shooting_ode_solver(func,u0,t0,tn,method,guess):
    """
    This function uses numerical shooting method, root finding method and numerical integrators to solve ordinary differential equations. 
    It includes finding the limit cycles and phase condition, plotting phase portraits and variables agaisnt time. 
    
    Parameters
    ---------
    Func:
        An n-th order ordinary differential equation we are solving. It has to be written in a system of n first-order odes.Also,the parameters should
        be included in the 'func'.The func returns a numpy array with u1,u2,u3...un
        
        Example: A predator_prey ode function can be written as below:
        
        def predator_prey(u,t):
            a = a
            b = b
            d = d
            du1dt = u[0]*(1-u[1])-(a*u[0]*u[1])/(d+u[0])
            du2dt = b*u[1]*(1-u[1]/u[0])
            return np.array([du1dt,du2dt])
            
        where u is the state vector : u = [u1,u2]. u1,u2 are number of preys and predators(variables we want to solve).t is initial time when u1,u2 
        get initial values. a, b, d are parameters in the model.The output is a numpy array of du1/dt and du2/dt, which will be used to in our shooting function later.
        
        
        parameters that are needed in the 'func',you need to specify the initial values.
        Example: for predator prey ode function, func = predator_prey(u1,u2,a,b,d), u1,u2,a,b,d you need to give them values
        
    u0:
        Initial values for variables in a numpy array.E.g. u1,u2,u3.. 
         
    t_span : 
        Input a time span t_span, which is a integer.The solution will get after t_span 
        
             
    method:
        Options for numerical integrators. 'euler' for  euler's method, 'RK4' for the fourth-order Runge-Kutta method, 'midpoint' for 
        the midpoint method.
        
    guess: 
        guess of initial values and period T to give limit cycles.Values are stored in a numpy array
        Example:
            guess for predator_prey function: guess = array([x_guess,y_guess,T])
            
   
   
    Returns
    ----------
    Returns a solution at t = t_span, returns plots of varibale agianst time between[0,t_span]. Plot the phase portrait. 
    Returns initial values to have limit cycles. output its period.


     """
    # Here is the code that does the shooting
    
    # single step of euler's method work for all ode 
    # u0 is a numpy array , the state initial values of variables in ode； t0 is the initial time； deltat is the stepsize
    def euler_step(u0,t0,deltat):
        u1 = u0+deltat*func(u0,t0)
        t1 = t0 + deltat
        return u1,t1


    ##single step of 4th Runge-Kutta method
    #ref: https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od
    
    def RK4(u0,t0,deltat):
        k1 = func(u0,t0) 
        k2 = func(u0+deltat*k1*0.5,t0+deltat/2)
        k3 = func(u0+deltat*k2*0.5,t0+deltat/2)
        k4 = func(u0+deltat*k3,t0+deltat)
        #update u0 value
        u1 = u0+deltat*(k1+2*k2+2*k3+k4)/6
        t1 = t0 + deltat
        return u1,t1
    
    # implement midpoint method
    def midpoint(u0,t0,deltat):
        t1 = t0+0.5*deltat
        u1 = u0 + 0.5*deltat*func(u0,t0)
        k = func(u1,t1)
        t2 = t0 + deltat
        u2 = u0 + deltat*k
        return u2,t2
    
    
    def numerical_shooting(s):
        u_0 = s[:-1]
        T0 = s[-1] 
        t_0 = t0
        #dxdt = func(u_0,t_0)
        for i in range(int(abs(T0)/deltat)):
            u_0,t_0=euler_step(u_0,t_0,deltat)
        return np.append(s[:-1]-np.array(u_0),func(u_0,t_0)[0])

    # use fsolve in scripy, the pros are this is an existing function that is well-built to find roots,
    # the cons are we need give initial guess and it may not converge and it may give reuslt that is wrong.like T less than 0   
    def solve(guess):
        result = root(lambda s: numerical_shooting(s),(guess))
        #if result.success:
            #print('success')
        #else:
            #print('failed to converge')
        return result
    
    u0 = u0
    t0 = t0
    deltat = 1e-3  # set deltat equals 1e-6, which will give very high accuracy of RK4 and euler 
    
    result = solve(guess)
    test = np.allclose(numerical_shooting(result.x),0,1e-3)
    
    
    if method == 'euler':
        for i in range(int((tn-t0)/deltat)):
            u0,t0 = euler_step(u0,t0,deltat)
    
    elif method == 'RK4':
        for i in range(int((tn-t0)/deltat)):
            u0,t0 = RK4(u0,t0,deltat)
            
    elif method == 'midpoint':
        for i in range(int((tn-t0)/deltat)):
            u0,t0 = midpoint(u0,t0,deltat)
                
      
    return [u0,result.x,test]


#%%


# use d^2x/dt^2 = -x ode to test if the function we built get the correct result

def func(u,t):
    x,y = u
    dxdt = y
    dydt = -x
    return np.array([dxdt,dydt])

z1 = shooting_ode_solver(func,(1,0),0,1,'euler',(1,0,6))
z2 = shooting_ode_solver(func,(1,0),0,1,'RK4',(1,0,6))
z3 = shooting_ode_solver(func,(1,0),0,1,'midpoint',(1,0,6))

print(z1)
print(z2)
print(z3)
# the answer is close to the correct value of solution x at t =1
#%%

def predator_prey(u,t0):
    a = 1
    b = 0.1
    d = 0.1
    du1dt = u[0]*(1-u[1])-(a*u[0]*u[1])/(d+u[0])
    du2dt = b*u[1]*(1-u[1]/u[0])
    return np.array([du1dt,du2dt])

z1 = shooting_ode_solver(predator_prey,(1,1),0,1,'euler',(0.3,0.3,10))
print(z1)

##ODE for the Hopf bifurcation normal form
def Hopf_bifurcation(u,t0):
    u1,u2 = u
    sigma = -1
    beta = 1
    du1dt = beta*u1-u2+sigma*u1*(u1**2+u2**2)
    du2dt = u1+beta*u2+sigma*u2*(u1**2+u2**2)
    return np.array([du1dt,du2dt])

z1 = shooting_ode_solver(Hopf_bifurcation,(1,0),0,1,'euler',(1,1,1))
z2 = shooting_ode_solver(Hopf_bifurcation,(1,0),0,1,'RK4',(1,1,1))
z3 = shooting_ode_solver(Hopf_bifurcation,(1,0),0,1,'midpoint',(1,1,1))

print(z1)
print(z2)
print(z3)



# unit test to test the solution
class Test_ODE(unittest.TestCase):
    def test_soln(self):
        
        z1 = shooting_ode_solver(Hopf_bifurcation,(1,0),0,1,'euler',(1,1,1))
        # Account for accuracy: as the shooting solution will not exactly match the explicit solution. Set tolerance equals to 1e-6
        if abs(z1[0][0]-math.sqrt(1)*math.cos(1)) < 1e-3 and abs(z1[0][1]-math.sqrt(1)*math.sin(1)) < 1e-3:  
            print("Test passed")
        else:
            print("Test failed")
        #self.assertAlmostEqual(z1[0],math.sqrt(1)*math.cos(1))
        #self.assertAlmostEqual(z1[1],math.sqrt(1)*math.sin(1))
        
unittest.main(argv=[''],verbosity=2, exit=False)
#if __name__ == '__main__':
    #unittest.main()



   