# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 18:39:32 2023

@author: YHU
"""

import numpy as np
import matplotlib.pyplot as plt

def predator_prey(x0,y0,a,b,d,deltat):
    dxdt = x0*(1-x0)-(a*x0*y0)/(d+x0)
    dydt = b*y0*(1-y0/x0)
    x1 = x0 + dxdt*deltat
    y1 = y0 + dydt*deltat
    return x1,y1


def predator_prey_plot(b):
        
    a = 1
    d = 0.1 
    deltat = 0.01
    x0,y0 = 0.05,0.05
    solution_x = [x0]
    solution_y = [y0]
    t = np.linspace(0,200,int(200/deltat)+1)
    for i in range(int(200/deltat)):
        x0,y0 = predator_prey(x0,y0,a,b,d,deltat)
        solution_x.append(x0)
        solution_y.append(y0)

    plt.figure()
    plt.plot(t,solution_x)
    plt.plot(t,solution_y)
    plt.legend(['Prey','Predator'])
    plt.yticks(np.arange(0,1,0.1))
    plt.xlabel('Time')
    plt.ylabel('Population')
    return plt




predator_prey_plot(0.4)
predator_prey_plot(0.2)


