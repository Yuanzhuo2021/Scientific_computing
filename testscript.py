# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 19:59:13 2023

@author: YHU
"""

import unittest

class Test_ODE(unittest.TestCase):
    def test_soln(self):
        self.assertAlmostEquals(1,1)
        
unittest.main(argv=[''],verbosity=0, exit=False)