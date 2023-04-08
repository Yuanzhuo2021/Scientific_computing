# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 23:53:31 2023

@author: YHU
"""
import numpy as np
import matplotlib.pyplot as plt

def shooting_ode_solver(func,u0,t_span,method,guess):
    """
    This function uses numerical shooting method and root finding method to solve ordinary differential equations. 
    It includes finding the limit cycles and phase condition, plotting phase portraits and variables agaisnt time. 
    
    Parameters
    ---------
    Func:
        An n-th order ordinary differential equation we are solving. It has to be written in a system of n first-order odes.Also,the parameters should
        be included in the 'func'.The func returns a numpy array with u1,u2,u3...un
        
        Example: A predator_prey ode function can be written as below:
        
        def predator_prey(u1,u2,a,b,d):
            du1dt = u1*(1-u1)-(a*u1*u2)/(d+u1)
            du2dt = b*u2*(1-u2/u1)
            return np.array([du1dt,du2dt])
            
        where u1,u2 are number of preys and predators(variables we want to solve), a, b, d are parameters in the model.The output is a numpy
        array of du1/dt and du2/dt, which will be used to in our shooting function later.
        
        
        parameters that are needed in the 'func',you need to specify the initial values.
        Example: for predator prey ode function, func = predator_prey(u1,u2,a,b,d), u1,u2,a,b,d you need to give them values
         
    t : 
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
    
    
    
    # single step of euler's method 
    def euler_step(x0,t0,deltat):
        x1 = x0+deltat*func
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
    
    
    if method == 'Euler':
        while i <= num_steps:
            x1 = euler_step(x1,t1,deltat)
            t1 = t1+deltat
            i += 1
        break
                    
    elif method == 'RK4':
        while i <= num_steps:
            x1 = RK4(x1,t1,deltat)
            t1 = t1+deltat
            i += 1
        break
                    
    else:
        print('Wrong input,please try again')


    
    
        
        
        
         
    
        
        
        
        