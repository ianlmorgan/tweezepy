import os
import numpy as np
import unittest
from tweezepy.allanvar import *

this_dir = os.path.dirname(os.path.abspath(__file__))

class test_m_generator(unittest.TestCase):
    def setUp(self):
        self.N = 128
    def test_all(self):
        maxn = self.N//2
        expected_all = np.linspace(1,maxn,maxn,dtype='int')
        m = m_generator(self.N,taus = 'all')
        np.testing.assert_allclose(m,expected_all)
    def test_octave(self):
        maxn = int(np.log2(self.N/2)) 
        expected_octave = np.logspace(0,maxn-1,maxn,base=2,dtype='int')
        m = m_generator(self.N,taus = 'octave')
        np.testing.assert_allclose(m,expected_octave)
    def test_decade(self):
        maxn = int(np.floor(np.log10(self.N/2)))
        expected_decade = [np.array([1,2,4])*k for k in np.logspace(0,maxn,maxn+1,base=10,dtype='int')]
        expected_decade = np.ravel(expected_decade)
        m = m_generator(self.N,taus = 'decade')
        np.testing.assert_allclose(m,expected_decade)

class test_avar(unittest.TestCase):
    def setUp(self):
        path = os.path.join(this_dir,'data/trace.txt')
        self.trace = np.loadtxt(path) 

    def test_avar_standard(self):
        avdata = avar(self.trace,100,overlapping = False)
        path = os.path.join(this_dir,'data/avdata_standard.txt')
        test_data = np.loadtxt(path)        
        np.testing.assert_allclose(avdata,test_data)
        
    def test_avar_overlapping(self):
        avdata = avar(self.trace,100,overlapping=True)
        path = os.path.join(this_dir,'data/avdata_overlapping.txt')
        test_data = np.loadtxt(path)
        np.testing.assert_allclose(avdata,test_data)

    def test_avar_total(self):
        avdata = totvar(self.trace,100)
        path = os.path.join(this_dir,'data/avdata_total.txt')
        test_data = np.loadtxt(path)
        np.testing.assert_allclose(avdata,test_data)

if __name__ == '__main__':
    unittest.main()