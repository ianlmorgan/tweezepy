# -*- coding: utf-8 -*-
"""
Created on Wed May 27 14:59:09 2020

Author: Ian L. Morgan
email: ilmorgan@ucsb.edu
"""
import math
import numpy as np

from tqdm.auto import tqdm

if __name__ == '__main__':
    
    F = 2           # force applied to bead in pN
    Z0 = 1000       # extension of molecule in nm
    
    alpha = 1.5e-5  # dissipation due to viscous drag, in pN s/nm; 
                    # 1.5e-5 is a typical value for an MT experiment 
                    # using a 1-micron bead
    sim_freq = 3e5  # frequency, in Hz, at which the simulation will be 
                    # conducted; this should be several orders of magnitude 
                    # faster than the downsampling frequency, below        
    freq = 400      # simulated data will be downsampled by averaging to this 
                    # frequency, akin to the camera frame rate   
    sim_time = 10   # length of simulation, in s
    
    num_sims = 1000 # number of simulations to be conducted; the results of 
                    # these will be combined to give estimates of uncertainty 
                    # on the fit parameters
                    
    kT = 4.1        # thermal energy, in pN nm (approx. 4.1 at room temp.)
    
    trace_length = int(sim_time*sim_freq) # number of points per simulation

    dt = 1/sim_freq    #length of time corresponding to each data point

    kx = F/Z0    # spring constant in x
                 # .002 pN/nm
    
    X0 = 0 # avg. x-position of bead
    
    bin_size = math.floor(sim_freq/freq)
    num_bins = math.floor(trace_length/bin_size)
    time_ds = np.arange(1,num_bins+1)/freq
    
    for n in tqdm(range(num_sims)):
        F_L = np.random.standard_normal(trace_length)
        xtrace = np.zeros(trace_length)
        for i in range(trace_length-1):
            dx =  dt/alpha*(math.sqrt(2*alpha*kT/dt)*F_L[i+1]-kx*(xtrace[i]-X0))
            xtrace[i+1] = xtrace[i] + dx
        xtrace_ds = np.zeros(num_bins)
        for i in range(num_bins):
            xtrace_ds[i] = np.mean(xtrace[(i*bin_size):(i+1)*bin_size])
        trace_ds = np.vstack((time_ds,xtrace_ds))
        np.savetxt("data/simulations/sim%s.csv"%n, trace_ds, delimiter=",")