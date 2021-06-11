import json
import numpy as np
import unittest

import tweezepy.cd

from matplotlib.testing.decorators import cleanup

class FunctionTests(unittest.TestCase):
    def test_loglikelihood(self):
        self.assertEqual(smm.loglikelihood(1,1,1),-1.0)

    def test_aliasPSD(self):
        self.assertAlmostEqual(smm.aliasPSD(1,1,1,1),8.872208996328476)
    
    def test_lansdorpPSD(self):
        self.assertAlmostEqual(smm.lansdorpPSD(1,1,1,1),8.2)
    
    def test_lansdorpPSD2(self):
        self.assertAlmostEqual(smm.lansdorpPSD2(1,1,1,1,1),9.2)

    def test_SMMAV(self):
        self.assertAlmostEqual(smm.SMMAV(1,1,1),1.3783481739415415)
    
    def test_SMMAV2(self):
        self.assertAlmostEqual(smm.SMMAV2(1,1,1,1),2.3783481739415415)
    
    def test_m_generator(self):
        # octave
        m = smm.m_generator(10,taus = 'octave')
        self.assertIsNone(np.testing.assert_array_equal(m, np.array([1,2])))
        # all
        m = smm.m_generator(10,taus = 'all')
        self.assertIsNone(np.testing.assert_array_equal(m, np.arange(1,11)))
        # decade
        m = smm.m_generator(10,taus = 'decade')
        self.assertIsNone(np.testing.assert_array_equal(m, np.array([1,2,4,10,20,40])))

class PSDTestCase(unittest.TestCase):
    def setUp(self):
        # Data simulated and downsampled to 100 Hz
        # gamma = 1.0E-5, kappa = 0.002
        npzdata = np.load('test-data/test_data.npz')
        self.xtrace = npzdata['data']
        self.psd = smm.PSD(self.xtrace,100,2)

    def test_psd(self):
        f,dens = smm.psd(self.xtrace,100,2)
        data = np.load('test-data/test_psd.npz')
        self.assertIsNone(np.testing.assert_array_equal(f, data['f']))
        self.assertIsNone(np.testing.assert_array_equal(dens, data['dens']))
    def test_PSD(self):
        data = self.psd.data
        expected_data = np.load('test-data/test_psdclass.npz')
        for key in data.keys():
            self.assertIsNone(np.testing.assert_array_equal(data[key], expected_data[key], key))

    def test_MLEfit(self):
        psd = self.psd
        psd.mlefit()
        results = psd.results
        expected = np.load('test-data/test_PSD_MLEfit.npz')
        for key in results.keys():
            self.assertIsNone(np.testing.assert_array_equal(results[key], expected[key]))
        # Incomplete testing
        psd.mlefit(fitfunc = 'lansdorpPSD',alpha = 1E-5)
        psd.mlefit(fitfunc = 'lansdorpPSD2')
        psd.mlefit(fitfunc = 'lansdorpPSD2',alpha = 1E-5)
        psd.mlefit(fitfunc = 'aliasPSD')
        psd.mlefit(fitfunc = 'aliasPSD', alpha = 1E-5)
    # Test figures
    @cleanup
    def test_create_figure(self):
        """
        very simple example test that creates a figure using pyplot.
        """
        psd = self.psd
        psd.mlefit()
        psd.plot()

    @cleanup
    def test_mcmc(self):
        psd = self.psd
        psd.mlefit()
        np.random.seed(0)
        psd.mcmc(walkers=10,steps =200)
        expected = np.load('test-data/test_PSD_mcmc.npz')
        self.assertIsNone(np.testing.assert_array_equal(psd.samples, expected['samples']))

class AVTestCase(unittest.TestCase):
    def setUp(self):
        # Data simulated and downsampled to 100 Hz
        # alpha = 1.0E-5, kappa = 0.002
        npzdata = np.load('test-data/test_data.npz')
        self.xtrace = npzdata['data']
        self.av = smm.AV(self.xtrace,100)

    def test_oavar(self):
        tau,eta,oav = smm.oavar(self.xtrace,100)
        data = {'tau':tau,'eta':eta,'oav':oav}
        expected = np.load('test-data/test_oavar.npz')
        for key in data.keys():
            self.assertIsNone(np.testing.assert_array_equal(data[key], expected[key]))
        
    def test_totvar(self):
        tau,eta,oav = smm.totvar(self.xtrace,100)
        data = {'tau':tau,'eta':eta,'oav':oav}
        expected = np.load('test-data/test_totvar.npz')
        for key in data.keys():
            self.assertIsNone(np.testing.assert_array_equal(data[key], expected[key]))

    def test_AV(self):
        data = self.av.data
        expected = np.load('test-data/test_AV.npz')
        for key in data.keys():
            self.assertIsNone(np.testing.assert_array_equal(data[key], expected[key]))

    def test_MLEfit(self):
        av = self.av
        av.mlefit()
        results = av.results
        expected = np.load('test-data/test_AV_MLEfit.npz')
        for key in results.keys():
            self.assertIsNone(np.testing.assert_array_equal(results[key], expected[key]))
        # Incomplete testing on fit methods
        av.mlefit(fitfunc='SMMAV',alpha = 1E-5)
        av.mlefit(fitfunc='SMMAV2')
        av.mlefit(fitfunc='SMMAV2',alpha = 1E-5)
    # Test figures
    @cleanup
    def test_create_figure(self):
        """
        very simple example test that creates a figure using pyplot.
        """
        av = self.av
        av.mlefit(pedantic = False)
        av.plot()

    @cleanup
    def test_mcmc(self):
        av = self.av
        av.mlefit()cd
        np.random.seed(0)
        av.mcmc(walkers=10,steps =200)
        expected = dict(np.load('test-data/test_AV_mcmc.npz'))
        self.assertIsNone(np.testing.assert_array_equal(av.samples, expected['samples']))

if __name__ == '__main__':
    unittest.main()