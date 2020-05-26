# -*- coding: utf-8 -*-
"""
Created on Fri May  8 09:49:47 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""

import pandas as pd

from pathlib import Path
from scipy.interpolate import UnivariateSpline

class Trace:
    def __init__(self,path,nrefs,nexps):
        """
        Reads txt file and creates a dataframe that can be modified.
        Assumes that the columns are in the labview data output format
        [Time,MP,MR,(Xref,Yref,Zref)*#ref,(Xexp,Yexp,Zexp)*#exp]
        """
        header = ['Time','MP','MR'] 
        header += ['%sref%s'%(p,bead) for bead in range(nrefs) for p in ('X','Y','Z')] 
        header += ['%sexp%s'%(p,bead) for bead in range(nexps) for p in ('X','Y','Z')]
        df = pd.read_csv(path,delimiter='\t',header = None)
        df.columns = header
        df = df[df.Time != 0]
        #df.update(df.Time.sub(df.Time.iloc[0]))
        #df = df.set_index('Time')
        
        self.df = df
    def correct(self,idxcorr,pxn):
        """
        Corrects x and y from pixels to nms. 
        Corrects z from uncorrected microns to nms.
        """
        dfcopy = self.df.copy() # create a copy of the dataframe
        xykeys = [x for x in dfcopy.keys() if 'X' in x or 'Y' in x]
        zkeys = [x for x in dfcopy.keys() if 'Z' in x]
        dfcopy[xykeys] = dfcopy[xykeys]*pxn
        dfcopy[zkeys] = dfcopy[zkeys]*idxcorr*1000
        self.df = dfcopy
        return self
    def dropbeads(self,refbeads,expbeads):
        """
        Creates new dataframe with only accepted reference and experimental beads.
        """
        header = ['Time','MP','MR'] 
        header += ['%sref%s'%(p,bead) for bead in refbeads for p in ('X','Y','Z')] 
        header += ['%sexp%s'%(p,bead) for bead in expbeads for p in ('X','Y','Z')]
        self.df = self.df[header]
        #self.xrkeys = ['Xref%s'%bead for bead in refbeads]
        #self.xekeys = ['Xexp%s'%bead for bead in expbeads]
        #self.yrkeys = ['Yref%s'%bead for bead in refbeads]
        #self.yekeys = ['Yexp%s'%bead for bead in expbeads]
        #self.zrkeys = ['Zref%s'%bead for bead in refbeads]
        #self.zekeys = ['Zexp%s'%bead for bead in expbeads]
        return self
    def refsubtraction(self,refbeads):
        """
        Subtracts the common-mode (average) of the reference beads in x, y, and z.
        """
        dfcopy = self.df.copy()
        for p in ['X','Y','Z']:
            rhead = [key for key in dfcopy.keys() if '%sref'%p in key]
            mean = dfcopy[rhead].mean(axis=1)
            keys = [key for key in dfcopy.keys() if p in key]
            dfcopy[keys] = dfcopy[keys].sub(mean,axis = 'index')
        self.df = dfcopy
        return self
    def process(self,refbeads,expbeads,idxcorr = 1.33/1.51,pxn=131.9):
        self.correct(idxcorr,pxn)
        self.dropbeads(refbeads,expbeads)
        self.refsubtraction(refbeads)
        return self
    def surface_correction(self,surfs):
        self.df.update(self.df.sub(surfs))
        return self
def find_surface(df):
    """
    Find the surface by choosing the 1% value at zero force.
    """
    surfs = df[[key for key in df.keys() if 'Zexp' in key]].quantile(0.05)
    surfs.name = 'Surface'
    return surfs
class Refidx:
    """
    Takes salt concentration in solution and returns refractive index correction.
    
    This functions creates and stores interpolation functions for the refractive
    index of various salt solutions. It assumes that you are using immersion oil
    with a refractive index of 1.515. 
    """
    def __init__(self):
        reffolder = Path('Refractive indices')
        refpaths = {path.with_suffix('').name:path for path in reffolder.iterdir()}

        spl = {}
        for key in refpaths.keys():
            refract = pd.read_csv(refpaths[key],header=None)
            refract.columns = ['C','n']
            spl[key] = UnivariateSpline(refract.C,refract.n)
        self.spl = spl
    def idxcorr(self,salt,conc):
        return self.spl[salt](conc/1000)/1.515