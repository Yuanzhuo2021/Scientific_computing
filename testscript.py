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
import continuation
import pde_solver
import finite_difference


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

def algebraic_cubic(u):
    x,c = u
    return x**3 - x + c


# test our codes
class TEST(unittest.TestCase):
    
 
    # test ode_solver
    def test_ode_solver(self):
        print('RK4 Pass: ',np.allclose(ode_solver.solve_to(func1,(1,1),0,1,'RK4',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1))))
        print('euler Pass: ',np.allclose(ode_solver.solve_to(func1,(1,1),0,1,'euler',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1)),1e-3))
        print('midpoint Pass: ',np.allclose(ode_solver.solve_to(func1,(1,1),0,1,'midpoint',0.001)[0],(math.sin(1)+math.cos(1),math.cos(1)-math.sin(1))))
    
    def test_shooting(self):
        z1 = shooting.shooting_solve(predator_prey,(1,1,1))
        print('shooting Pass: ', np.allclose(shooting.shooting_func(predator_prey,z1),[0,0,0],1e-5))

    def test_contiuation(self):
        z2 = continuation.natural(algebraic_cubic,(-2,-1.5),2)
        z3 = continuation.psuedo_arclength(algebraic_cubic,(-2,-1.5),2)
        print('conitunation Pass: ',np.allclose(z2[0],z3[0][0:-2],1e-3))
        
    def test_finite_difference(self):
        q = lambda x: x*0 +1
        z4 = finite_difference.fd(0,10,1,2,2,q)
        soln = (-1/(2*2))*(z4[1]-0)*(z4[1]-10)+((2-1)/(10-0))*(z4[1]-0)+1
        print('finite_difference Pass: ', np.allclose(z4[0],soln,1e-3))
        
    def test_pde_solver(self):
        f = lambda x: np.sin(math.pi*(x-0)/(1-0))
        z5 = pde_solver.explicit_euler(0,1,0,0,f,0.5,1,'euler')    
        real5 = np.exp(-0.5*((math.pi)**2)*(1)/((1-0)**2))*np.sin(math.pi*(z5[0]-0)/(1-0))    
        print('explicit_euler Pass: ',np.allclose(real5,z5[1][-1],1e-2))

        z6 = pde_solver.explicit_euler(0,1,0,0,f,0.5,1,'RK4')      
        real6 = np.exp(-0.5*((math.pi)**2)*(1)/((1-0)**2))*np.sin(math.pi*(z6[0]-0)/(1-0))   
        print('explicit_euler(RK4) Pass: ',np.allclose(real6,z5[1][-1],1e-2))

        z7 = pde_solver.implicit_euler(0,1,0,0,f,0.5,1)
        real7 = np.exp(-0.5*((math.pi)**2)*(1)/((1-0)**2))*np.sin(math.pi*(z7[0]-0)/(1-0))  
        print('implicit_euler Pass: ',np.allclose(real7,z7[1][-1],1e-1))
        
        z8 = pde_solver.crank_nicolson(0,1,0,0,f,0.5,1)
        real8 = np.exp(-0.5*((math.pi)**2)*(1)/((1-0)**2))*np.sin(math.pi*(z8[0]-0)/(1-0))  
        print('crank_nicolson Pass: ',np.allclose(real8,z8[1][-1],1e-4))
        
        
        
        
        
        
if __name__ == '__main__':
    #unittest.main()
    unittest.main(argv=[''],verbosity=2, exit=False)
    
    
    
    
    
    

    