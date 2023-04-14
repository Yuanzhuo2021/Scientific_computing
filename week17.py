# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 04:51:56 2023

@author: YHU
"""

import numpy as np
from scipy.optimize import fsolve,root
import matplotlib.pyplot as plt

def algebraic_cubic(x,c):
    return x**2 - x + c

def continuation(func,rang):
    S = 1000
    step = (rang[1]-rang[0])/S
    # first generate the first value of solution
    result= root(lambda u: func(u,rang[0]),-2)
    print(rang[0])
    print(result.x)
    solution = result.x
    C = [rang[0]]
    for i in range(S):
        c = rang[0] + step*(i+1)
        result = root(lambda u: func(u,c),result.x)
        if result.success:
            result = result
            C.append(c)
            solution = np.hstack((solution,result.x))
        else:
            print(result.message)
            print('Converge fails at c = ',c)
            break
        
    return [C,solution]

z = continuation(algebraic_cubic,(-2,2))



plt.plot(z[0],z[1])
plt.xlabel('Parameter')
plt.ylabel('Solution')