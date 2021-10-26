import os
import numpy as np
import unittest

from tweezepy import AV,PSD

from matplotlib.testing.decorators import cleanup
this_dir = os.path.dirname(os.path.abspath(__file__))

class test_calibration(object):
    def setUp(self):
        test_trace = os.path.join(this_dir,'data/trace.txt')
        self.trace = np.loadtxt(test_trace)
        self.fsample = 100
        self.method = self.method(self.trace,self.fsample)
        self.method.mlefit()
        self.test_trace()
        self.test_fsample()
    
    def test_trace(self):
        self.assertEqual(self.trace.size,self.method.trace.size)

    def test_fsample(self):
        self.assertEqual(self.fsample,self.method.fsample)

    def test_params(self):
        np.testing.assert_allclose(self.params.size,self.method.params.size)
    
    def test_std_errors(self):
        np.testing.assert_allclose(self.std_errors.size,self.method.std_errors.size)
    
    def test_results(self):
        self.test_params()
        self.test_std_errors()

    @cleanup
    def test_plot(self):
        self.method.plot()
    
class testPSD(unittest.TestCase,test_calibration):
    def setUp(self):
        self.method = PSD
        test_calibration.setUp(self)
        
        path = os.path.join(this_dir,'data/psd_params.txt')
        self.params = np.loadtxt(path)
        path = os.path.join(this_dir,'data/psd_std_errors.txt')
        self.std_errors = np.loadtxt(path)

    def test_psd_results(self):
        self.test_results()

    def test_psd_plot(self):
        self.test_plot()


class testAV(unittest.TestCase,test_calibration):
    def setUp(self):
        self.method = AV
        test_calibration.setUp(self)
        path = os.path.join(this_dir,'data/av_params.txt')
        self.params = np.loadtxt(path)
        path = os.path.join(this_dir,'data/av_std_errors.txt')
        self.std_errors = np.loadtxt(path)
        
    def test_av_results(self):
        self.test_results()

    def test_avplot(self):
        self.test_plot()

if __name__ == '__main__':
    unittest.main()