# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 19:59:13 2023

@author: YHU
"""

import unittest
import math
import numpy as np
import ode_solver
import shooting
#import continuation
#import implicit_euler
#import explicit_euler



def func1(u,t):
    x,y = u
    dxdt = y
    dydt = -x
    return np.array([dxdt,dydt])   

def predator_prey(u,t0):
    x,y = u
    a = 1
    b = 0.1
    d = 0.1
    dxdt = x*(1-x)-(a*x*y)/(d+x)
    dydt = b*y*(1-y/x)
    return np.array([dxdt,dydt])

# test our codes
class TEST(unittest.TestCase):
    
 
    # test ode_solver
    def test_ode_solver(self):
        np.testing.assert_array_almost_equal(ode_solver.solve_to(func1,(1,1),0,1,'RK4',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1)))
        np.testing.assert_array_almost_equal(ode_solver.solve_to(func1,(1,1),0,1,'euler',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1)),decimal=4)
        np.testing.assert_array_almost_equal(ode_solver.solve_to(func1,(1,1),0,1,'midpoint',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1)))
    
    def test_shooting(self):
        z1 = shooting.shooting_solve(predator_prey,(1,1,1))
        np.testing.assert_array_almost_equal(shooting.shooting_func(predator_prey,z1),[0,0,0],1e-5)

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    #unittest.main()
    unittest.main(argv=[''],verbosity=4, exit=False)
    
    
    
    
    
    

    