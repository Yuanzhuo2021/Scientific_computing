# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 04:51:56 2023

@author: YHU
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def algebraic_cubic(x,c):
    return x**3 - x + c

def continuation(func,rang):
    S = 1000
    step = (rang[1]-rang[0])/S
    # first generate the first value of solution
    result = fsolve(lambda u: func(u,rang[0]),2)
    print(rang[0])
    print(result)
    solution = result
    C = [rang[0]]
    for i in range(S):
        c = rang[0] + step*(i)
        C.append(c)
        result = fsolve(lambda u: func(u,c),result)
        solution = np.hstack((solution,result))
    return np.array([C,solution])

z = continuation(algebraic_cubic,(-2,2))


plt.plot(z[0],z[1])