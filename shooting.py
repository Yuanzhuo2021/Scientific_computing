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
import ode_solver

def shooting_func(func,s):
    """
    This function uses numerical shooting method and numerical integrators to represent shooting equations

    Parameters
    ---------
    Func:
        An n-th order ordinary differential equation we are solving. It has to be written in a system of n first-order odes.Also,the parameters should
        be included in the 'func'.The func returns a numpy array with u1,u2,u3...un
        
        Example: A predator_prey ode function can be written as below:
        
        def predator_prey(u,t):
            x,y = u
            a = a
            b = b
            d = d
            dxdt = x*(1-y)-(a*x*y)/(d+x)
            dydt = b*y*(1-y/x)
            return np.array([dxdt,dydt])
            
        where u is the state vector : u = (x,y). x,y are number of preys and predators(variables we want to solve).t is initial time when x,y
        get initial values. a, b, d are parameters in the model.The output is a numpy array of dx/dt and dy/dt, which will be used in our shooting function later.
        
        
    s: 
       The state vector and  period T to give limit cycles.Values are stored in a numpy array
        Example:
            guess for predator_prey function: guess : (x_guess,y_guess,T)
            
    Returns
    ----------
    list
    Returns [x(0)-x(T), y(0)-]


     """
    # Here is the code that does the shooting
    
    method = 'RK4'  # use RK4 because it has the highest accuracy among our methods
    try:
        # use solve_to function in ode_solver library. The solution contains solution at t = tn
        z = ode_solver.solve_to(func,s[:-1],0,s[-1],method,0.001)

        #return [s[0]-z[0][0],s[1]-z[0][1],func(s[:-1],0)[0]] # z[0][0] is x, z[0][1] is y
        return [s[0]-z[0][0],s[1]-z[0][1],func(s[:-1],0)[0]]
    except:
        print('Wrong input')
        return [] 
    
     
def shooting_solve(func,guess):
    """
    This function uses root finding method to solve shooting equation
    
    Parameters
    ---------
    Func: callable
        An n-th order ordinary differential equation we are solving. It has to be written in a system of n first-order odes.Also,the parameters should
        be included in the 'func'.The func returns a numpy array with u1,u2,u3...un
        
        Example: A predator_prey ode function can be written as below:
        
        def predator_prey(u,t):
            x,y = u
            a = a
            b = b
            d = d
            dxdt = x*(1-y)-(a*x*y)/(d+x)
            dydt = b*y*(1-y/x)
            return np.array([dxdt,dydt])
            
        where u is the state vector : u = (x,y). x,y are number of preys and predators(variables we want to solve).t is initial time when x,y 
        get initial values. a, b, d are parameters in the model.The output is a numpy array of dx/dt and dy/dt, which will be used to in our shooting function later.
        
    
    guess: numpy.array
        guess of initial values of varibles and period T to give limit cycles.Values are stored in a numpy array
        Example:
            guess for predator_prey function: guess : (x_guess,y_guess,T)
            
    Return
    -------
    numpy array. The state vector that will have limit cycles and its period
            

    """
     # use root in scipy, the pros are this is an existing function that is well-built to find roots,
    # the cons are we need give initial guess and it may not converge and it may give reuslt that is wrong.like T less than 0 
    result = root(lambda s: shooting_func(func,s),guess)
    if result.success:
        result = result.x
        
    else:
        result = np.nan
        print('failed to converge')
    return result


#%%   
if __name__ == "__main__":

    def predator_prey(u,t0):
        x,y = u
        a = 1
        b = 0.1
        d = 0.1
        dxdt = x*(1-x)-(a*x*y)/(d+x)
        dydt = b*y*(1-y/x)
        return np.array([dxdt,dydt])

    z1 = shooting_solve(predator_prey,(1,1,1))
    test = np.allclose(shooting_func(predator_prey,z1),0,1e-5) # test the result
    print('The shooting solution of predator prey is',z1)
    print(test)


#%%
    ##ODE for the Hopf bifurcation normal form
    def Hopf_bifurcation(u,t0):
        u1,u2 = u
        sigma = -1
        beta = 1
        du1dt = beta*u1-u2+sigma*u1*(u1**2+u2**2)
        du2dt = u1+beta*u2+sigma*u2*(u1**2+u2**2)
        return np.array([du1dt,du2dt])

    z2 = shooting_solve(Hopf_bifurcation,(1,1,1))
    test = np.allclose(shooting_func(Hopf_bifurcation,z2),0,1e-5)
    print('The shooting solution of Hopf bifurcation is',z2)
    print(test)




   