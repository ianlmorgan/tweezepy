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
    
    def test_trace(self):
        self.assertEqual(self.trace.size,self.method.trace.size)

    def test_fsample(self):
        self.assertEqual(self.fsample,self.method.fsample)

    def test_x(self):
        self.assertEqual(self.x.size,self.method.data['x'].size)

    def test_shape(self):
        self.assertEqual(self.shape.size,self.method.data['shape'].size)

    def test_y(self):
        self.assertEqual(self.y.size,self.method.data['y'].size)
    
    def test_data(self):
        self.test_x()
        self.test_shape()
        self.test_y()

    def test_params(self):
        self.assertEqual(self.params.size,self.method.params.size)
    
    def test_std_errors(self):
        self.assertEqual(self.std_errors.size,self.method.std_errors.size)
    
    def test_residuals(self):
        self.assertEqual(self.residuals.size,self.method.residuals.size)

    def test_results(self):
        self.test_params()
        self.test_std_errors()
        self.test_residuals()

    @cleanup
    def test_plot(self):
        self.method.plot()
    
class testPSD(unittest.TestCase,test_calibration):
    def setUp(self):
        test_calibration.setUp(self)
        self.method = PSD(self.trace,self.fsample,bins = 15)
        self.method.mlefit()
        self.test_trace()
        self.test_fsample()
        path = os.path.join(this_dir,'data/psd_f.txt')
        self.x = np.loadtxt(path)
        path = os.path.join(this_dir,'data/psd_shapes.txt')
        self.shape = np.loadtxt(path)
        path = os.path.join(this_dir,'data/psd_dens.txt')
        self.y = np.loadtxt(path)
        path = os.path.join(this_dir,'data/psd_params.txt')
        self.params = np.loadtxt(path)
        path = os.path.join(this_dir,'data/psd_std_errors.txt')
        self.std_errors = np.loadtxt(path)
        path = os.path.join(this_dir,'data/psd_residuals.txt')
        self.residuals = np.loadtxt(path)

    def test_psd_data(self):
        self.test_data()
    
    def test_psd_results(self):
        self.test_results()

    def test_psd_plot(self):
        self.test_plot()


class testAV(unittest.TestCase,test_calibration):
    def setUp(self):
        test_calibration.setUp(self)
        self.method = AV(self.trace,self.fsample)
        self.method.mlefit()
        self.test_trace()
        self.test_fsample()
        path = os.path.join(this_dir,'data/av_taus.txt')
        self.x = np.loadtxt(path)
        path = os.path.join(this_dir,'data/av_shapes.txt')
        self.shape = np.loadtxt(path)
        path = os.path.join(this_dir,'data/av_oavs.txt')
        self.y = np.loadtxt(path)
        path = os.path.join(this_dir,'data/av_params.txt')
        self.params = np.loadtxt(path)
        path = os.path.join(this_dir,'data/av_std_errors.txt')
        self.std_errors = np.loadtxt(path)
        path = os.path.join(this_dir,'data/av_residuals.txt')
        self.residuals = np.loadtxt(path)

    def test_avdata(self):
        self.test_data()
        
    def test_av_results(self):
        self.test_results()

    def test_avplot(self):
        self.test_plot()


if __name__ == '__main__':
    unittest.main()