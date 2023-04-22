# -*- coding: utf-8 -*-
"""
Implicit Euler

@author: YHU
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def implicit_euler(a,b,alpha,beta,D,t):
    """
    This function uses implicit euler method to solve pde

    Parameters
    ----------
    a : int
        Left boundary 
    b : int
        Right boundary 
    alpha : float
        Left boundary condition e.g. u(a,0) = alpha 
    beta : float
        Right boundary condition e.g. u(b,0) = beta
    D : float/int
        Parameter in pde
    t : int
        Final time solve to

    Returns
    -------
    numpy.array
    The solution u(x,t) of pde

    """
    
    # set number of steps of space and time
    N_space = 200
    N_time = 1000
    
    #step size
    dx = (b-a)/N_space
    dt = (t-0)/N_time
    
    C = D*dt/(dx**2)

    # identity matrix
    I = np.identity(N_space-1)
    
    #stepwise x values
    x = np.linspace(a,b,N_space+1)

    # use numpy to solve the system of equations. There are N-1 equations as there are N-1 unknowns.

    # build A-DD
    A = np.zeros((N_space-1,N_space-1))
    # the first row and the last row are quite different.We deal with them seperately
    A[0,0]=-2
    A[0,1]=1
    A[N_space-2,N_space-2]=-2
    A[N_space-2,N_space-3]=1
    # for the rest of the row, coefficient for u[i+1] is 1, u[i] is -2,u[i-1] is 1
    for i in range(N_space-3):
        A[i+1,i]=1
        A[i+1,i+1]=-2
        A[i+1,i+2]=1
    print(A)
    
    #build b_DD
    b_DD = np.zeros(N_space-1)
    #put initail and end conditions
    b_DD[0]=alpha
    b_DD[N_space-2]=beta
    print(len(b_DD))


    # implement the equation

    print(I - C*A)
    print(b_DD)

    U = []
    u0 = np.zeros(N_space-1)
    for n in range(N_time):
        u1 = np.linalg.solve(I-C*A,u0+C*b_DD)
        U.append(u1)
        u0 = u1 
        
    # add boundary conditions
    U = np.hstack((U,beta*np.ones((N_time,1)))) 
    U = np.hstack((alpha*np.ones((N_time,1)),U))
    print(U)
    
    return [x,U]


z = implicit_euler(0,1,1,0,0.5,1)


plt.plot(z[0],z[1][0])
plt.plot(z[0],z[1][200])
plt.plot(z[0],z[1][400])
plt.plot(z[0],z[1][999])
plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.legend(['t=0','t=0.2','t=0.4','t=1'])
