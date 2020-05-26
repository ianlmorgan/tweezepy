# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:53:19 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""

import pandas as pd

from tweepy.smmcalibration import AV

def fext_inner(f,plot = False):
    """Reads in csv file and produces dataframe with relevant fext values."""
    df = pd.read_csv(f)
    expbeads = [int(col.replace('Zexp','')) for col in df.columns if 'Zexp' in col]
    results = {}
    for bead in expbeads:
        results[bead] = {}
        
        xtrace = df['Xexp%s'%bead].to_numpy()
        xav = AV(xtrace,400)
        alpha,kappa = xav.mlefit()
        if plot == True:
            xav.plot()
        
        zmean = df['Zexp%s'%bead].mean() # mean length in nm
        force = zmean*kappa # force in pN
        
        results[bead] = {}
        results[bead]['force_pN'] = force
        results[bead]['length_nm_mean'] = zmean
        results[bead]['length_nm_std'] = zmean
        results[bead]['mp_mm'] = df['MP'].median()
        results[bead]['alpha'] = alpha
        results[bead]['kappa'] = kappa
    results = pd.DataFrame.from_dict(results,orient='index')
    return results

def fext(folder):
    """Reads in folder of trace files in csv format and produces dataframe with all fext values."""
    dfs = [fext_inner(f) for f in folder.glob('*.csv')]
    if not len(dfs) == 0:
        result = pd.concat(dfs,keys=range(len(dfs)),names = ['file','bead'])
        force_mean = result.groupby(['bead','mp_mm']).mean()['force_pN']
        result = result.join(force_mean,on=['bead','mp_mm'],rsuffix = '_mean')
        result = result.swaplevel(0,1)
        #result.plot('force_pN_avg','z_nm_mean',marker='o',logx=True,legend = False)
        #plt.show()
        return result