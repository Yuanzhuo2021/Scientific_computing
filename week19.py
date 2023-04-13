# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 21:59:11 2023

@author: YHU
"""

## finite differenting

import numpy as np


# convert ode BVP into a system of algebraic equations
a = 0 
b = 10
gamma1 = 1
gamma2 = 2
N =1000

dx = (b-a)/N


u_a = gamma1
u_b = gamma2

# second order differentiation: (u[i+1] -2u[i] + u[i-1])/(dx)^2 +q[i]= 0 .Because the source equation is 0, q[i] = 0.

A = np.zeros((N+1,N+1))

print(A)