# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 19:59:13 2023

@author: YHU
"""

import unittest
import math
import numpy as np
import ode_solver
#import shooting
#import continuation
#import implicit_euler
#import explicit_euler



def func2(u,t):
    x,y = u
    dxdt = y
    dydt = -x
    return np.array([dxdt,dydt])   

# test our codes
class TEST(unittest.TestCase):
    
 
    # test ode_solver
    def test_ode_solver(self):
        np.testing.assert_array_almost_equal(ode_solver.solve_to(func2,(1,1),0,1,'RK4',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1)))
        np.testing.assert_array_almost_equal(ode_solver.solve_to(func2,(1,1),0,1,'euler',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1)),decimal=4)
        np.testing.assert_array_almost_equal(ode_solver.solve_to(func2,(1,1),0,1,'midpoint',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1)))
    
    def test_shooting(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)

if __name__ == '__main__':
    #unittest.main()
    unittest.main(argv=[''],verbosity=4, exit=False)
    
    
    
    
    
    
  #%% 
    
    
    
    
    # unit test to test the solution
class Test_ODE(unittest.TestCase):
    def test_soln(self):
        
        z1 = shooting(Hopf_bifurcation,(1,0),0,1,'euler',(1,1,1))
        # Account for accuracy: as the shooting solution will not exactly match the explicit solution. Set tolerance equals to 1e-6
        if abs(z1[0][0]-math.sqrt(1)*math.cos(1)) < 1e-3 and abs(z1[0][1]-math.sqrt(1)*math.sin(1)) < 1e-3:  
            print("Test passed")
        else:
            print("Test failed")
        #self.assertAlmostEqual(z1[0],math.sqrt(1)*math.cos(1))
        #self.assertAlmostEqual(z1[1],math.sqrt(1)*math.sin(1))
        
unittest.main(argv=[''],verbosity=2, exit=False)
#if __name__ == '__main__':
    #unittest.main()
    