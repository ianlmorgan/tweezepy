"""
This file tests the simulations module.
"""
import numpy as np
import os
import unittest

from tweezepy.simulations import downsampled_trace

this_dir = os.path.dirname(os.path.abspath(__file__))
class test_simulations(unittest.TestCase):
    def setUp(self):
        path = os.path.join(this_dir,'data/trace.txt')
        self.test_trace = np.loadtxt(path)
    def test_downsampled_trace(self):
        test_trace = downsampled_trace(seed=0)
        np.testing.assert_allclose(self.test_trace,test_trace)

if __name__ == '__main__':
    unittest.main()
