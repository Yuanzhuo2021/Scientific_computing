# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 21:59:11 2023

@author: YHU
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as ss

def finite_difference(a,b,alpha,beta,D,q):
    """
    This function is a BVP solver that is capable of finding 
    numerical solutions to ODEs using finite difference method.

    Parameters
    ----------
    a : int/float
        Left boundary
    b : int/float
        Right boundary
    alpha : int/float
        Value at Left boundary 
    beta : int/float
        Value at Left boundary 
    D: int/float
        Parameter in the ODE
    q : callable
        User defined function of q(x)
        e.g. q = lambda x: 2*x+1 or 
            q = lambda x : x*0 + p (if q(x) = p for all x)

    Returns
    -------
    u : numpy array
        The solutions of ODE in the domain

    """
    
    # convert ode BVP into a system of algebraic equations

    N = 100000
    dx = (b-a)/N
    x = np.linspace(a,b,N-1)
    # second order differentiation: D*(u[i+1] -2u[i] + u[i-1])/(dx)^2 +q[i]= 0 

    # Build matrix vector form to represent the system of equations. There are N-1 equations as there are N-1 unknowns.
    # build A-DD
    A = ss.lil_matrix((N-1, N-1))

    # the first row and the last row are quite different.We deal with them seperately
    A[0,0]=-2
    A[0,1]=1
    A[N-2,N-2]=-2
    A[N-2,N-3]=1

    # for the rest of the row, coefficient for u[i+1] is 1, u[i] is -2,u[i-1] is 1
    for i in range(N-3):
       A[i+1,i]=1
       A[i+1,i+1]=-2
       A[i+1,i+2]=1

    #right side of the equations:b-DD
    B = np.zeros(N-1)
    #put initail and end conditions
    B[0]=alpha
    B[N-2]=beta


    # account for q term
    Q  = q(x)

    # use sparse matrix to solve the matrix-vector form of equation
    u = ss.linalg.spsolve(A, -B-(((dx)**2)*Q)/D)

    return [u,x]
# %%

if __name__ == '__main__':

    a = 0 
    b = 10
    alpha = 1
    beta = 2
    N =100000
    D =2
    q = lambda x: x*0 +1 # q[x] =1
    
    z = finite_difference(a,b,alpha,beta,D,q)
    
    # test the accuray, use finite difference, it equals real solution with an error tolerance 1e-3
    def real_soln2(x):
        soln = (-1/(2*D))*(x-a)*(x-b)+((beta-alpha)/(b-a))*(x-a)+alpha
        return soln
    print(np.allclose(z[0],real_soln2(z[1]),1e-3))

    plt.plot(z[1],real_soln2(z[1]),z[1],z[0])
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend(['real solution','Finite difference'])

