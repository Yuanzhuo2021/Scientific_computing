# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 04:51:56 2023

@author: YHU

continuations
"""

import numpy as np
from scipy.optimize import fsolve,root
import matplotlib.pyplot as plt


def natural(func,rang,guess):
    """
    This function is used natural parameter continuation method, get the solution depending on 
    parameter c

    Parameters
    ----------
    func : callable
        Equation to solve
        
        E.g.
        def algebraic_cubic(u):
            x,c = u
            return x**3 - x + c
        
    rang : numpy.array
        the range of parameter c
    guess : int/numpy.array
        Initial guess of the solution

    Returns
    -------
    Return a list containing values of solutions and corresponding parameters
        

    """
    S = 1000
    step = (rang[1]-rang[0])/S
    C = []
    solution = []
    for i in range(S):
        # simply add a stepsize to c 
        c = rang[0] + step*(i)
        result = root(lambda x: func((x,c)),guess)
        if result.success:
            C.append(c)
            solution = np.hstack((solution,result.x))
            result = guess
        else:
            print(result.message)
            print('Converge fails at c = ',c)
            break
    return [solution,C]




def psuedo_arclength(func,rang,guess):
    """
    This function is used to find solutions of an equation or systems of equations 
    that depend on a parameter c

    Parameters
    ----------
    func : callable
        ode function or equation
        
        Example:
        def algebraic_cubic(u):
            x,c = u
            return x**3 - x + c
        
    rang : numpy.array
        the range of parameter c
    guess : numpy.array/int
        Intial guess of the solution

    Returns
    -------
    Return a list with the values of the parameter and corresponding solutions

    """
    
    # use natural parameter continuation generate first two solutions 
    z = natural(func,rang,guess)
    if len(z[0])>=2 and len(z[1])>=2:
        # get the first two solns using natural parameter continuation
        u0 = np.array([z[0][0],z[1][0]])
        u1 = np.array([z[0][1],z[1][1]])  # solution, C
        
        
       # print(delta) 
        soln = [z[0][0],z[0][1]]
        param = [z[1][0],z[1][1]]
        
        # define an extra pseudo arclenth equation which is need to be solved
        def psuedo_equation(u,u2):
            return [np.dot((u-u2),delta),func(u)]
        
        for i in range(1000):
            delta = (u1-u0)
            u2 = u1 + delta # estimated value
            true_value = root(lambda u: psuedo_equation(u,u2),u2) # true value of u2

            if true_value.success:
                soln.append(true_value.x[0])
                param.append(true_value.x[-1])
                u0 = u1
                u1 = true_value.x
                
            else:
                print('converge fails at c =',u1[-1])
                break
    else:
        print('Converges fail, please choose a new guess')
    return [soln,param]

#%%

if __name__ == '__main__':
    
    # define a equation
    def algebraic_cubic(u):
        x,c = u
        return x**3 - x + c
    
    #solve equation using natural continuation
    z = natural(algebraic_cubic,(-2,2),2)
    plt.plot(z[1],z[0])
    plt.xlabel('Parameter')
    plt.ylabel('Solution')
    
    # solve using pseudo-arclength continuation
    zz = psuedo_arclength(algebraic_cubic,(-2,2),2)
    #print(zz[0])
    plt.figure()
    plt.plot(zz[1],zz[0])
    plt.xlabel('Parameter')
    plt.ylabel('Solution')

#%%





