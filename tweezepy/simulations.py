"""
This module include functions for running bead trajectory simulations.

"""
import numpy as np
from numba import jit

@jit(nopython=True) # speeds up certain types of code
def simulate_trace(gamma = 1.e-5,
                   kappa = None,
                   fsim = 400,
                   sim_points = 10240,
                   seed = None):
    """Takes parameters and gives a simulated bead trajectory.

    .. math::
        \Delta x = \\frac{\Delta t}{\\alpha}\left(\sqrt{\\frac{2\\alpha k_BT}{dt}} F_L-\\kappa x_{i-1}\\right)
    .. math:: 
        x_i = x_{i-1}+dx

    Parameters
    ----------
    gamma : float, optional
        Stokes dissipation in pNs/nm^2, by default 1.e-5
    kappa : float, optional
        Spring constant in pN/nm, by default None
    fsim : int, optional
        Simulation frequency in Hz, by default 100
    sim_points : int, optional
        Number of bead positions, by default 10240
    seed : int, optional
        Random seed, by default None
    
    Returns
    -------
    xtrace : array
        Array of bead positions.

    """
    if seed != None:
        np.random.seed(seed)
    if not kappa:
        kappa = 2*np.pi*10*gamma
    if fsim < kappa/gamma: # Simulation is unstable below this condition
        raise ValueError(r'The simulation frequency must be larger than kappa/gamma.')
    dt = 1/fsim # time step between collecting probe position
    kT = 4.1    # thermal energy, in pN nm (approx. 4.1 at room temp.)
    std = np.sqrt(2*gamma*kT/dt) # standard deviation of random gaussian force
    F_L = np.random.normal(0,std,sim_points-1) # Random Langevin force
    xtrace = np.zeros(sim_points) 
    for i in range(1,sim_points):
        dx =  dt/gamma * (F_L[i-1]-kappa*xtrace[i-1]) # Smoluchowski eqn
        xtrace[i] = xtrace[i-1] + dx
    return xtrace

def downsampled_trace(gamma = 1.e-5,
                      kappa = None,
                      fsample = 100,
                      N = 10240, 
                      seed = None):
    """
    Simulates and downsamples trace.

    Parameters
    ----------
    gamma : float, optional
        Stokes dissipation in pNs/nm^2, by default 1.e-5
    kappa : float, optional
        Spring constant in pN/nm, by default None
    fsample : int, optional
        Sampling frequency in Hz, by default 100
    N : int, optional
        Number of downsampled bead positions, by default 10240
    seed : int, optional
        Random seed, by default None

    Returns
    -------
    xtrace_ds : numpy.array
        Array of bead positions.

    """
    bin_size = 1000 # points to average when downsampling
    fsim = int(bin_size*fsample) # simulation frequency
    sim_points = N*bin_size # number of points to simulate
    xtrace = simulate_trace(gamma,kappa,fsim,sim_points,seed) # do simulation
    bins = xtrace.reshape((N,bin_size)) # reshape into bins 
    xtrace_ds = np.mean(bins,axis=1) # average bins to produce downsampled trace
    return xtrace_ds