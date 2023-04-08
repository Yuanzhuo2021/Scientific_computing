# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 23:53:31 2023

@author: YHU
"""

def shooting_ode_solver(func,u0,t_span,method,guess):
    """
    This function uses numerical shooting method and root finding method to solve ordinary differential equations. 
    It includes finding the limit cycles and phase condition, plotting phase portraits and variables agaisnt time. 
    
    Parameters
    ---------
    Func:
        An ordinary differential equation we are solving. It has to be written in a system of first order of odes.
        
        Example: A predator_prey ode function can be written as below:
        
        def predator_prey(x0,y0,a,b,d):
            dxdt = x0*(1-x0)-(a*x0*y0)/(d+x0)
            dydt = b*y0*(1-y0/x0)
            return np.array([dxdt,dydt])
            
        where x0,y0 are number of preys and predators(variables we want to solve), a, b, d are parameters in the model.The output is a numpy
        array of dx/dt and dy/dt, which will be used to in our shooting function later.
        
        
    u0:
        parameters that are needed in the ode function,stored in a numpy array.You need to specify the values of elements in u0.
        Example: for predator prey ode function, u0=array([x0,y0,a,b,d])
        
    t_span:
        
        
    method:
        Options for numerical integrators. 'euler' for  euler's method, 'RK4' for the fourth-order Runge-Kutta method, 'midpoint' for 
        the midpoint method.
        
    guess: 
        guess of initial values and period T to give limit cycles.Values are stored in a numpy array
        Example:
            guess for predator_prey function: guess = array([x_guess,y_guess,T])
            
   
   
    Returns
    ----------
    Return
    
    