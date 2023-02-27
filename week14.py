# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#use euler's method to solve a ode probelm
import numpy as np
import matplotlib as plt
import math



def euler_step(x0,y0,deltat):
    y = y0+deltat*x0
    return y



Y = exp(x)
x0 =1
deltat = 0.1
num_steps = 1/deltat
while 1 <= 1/deltat:
    y = euler_step(x0,deltat)
    x0 = x+deltat
    error = Y - y
    

print(y)
