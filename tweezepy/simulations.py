import numpy as np
from numba import jit
@jit(nopython=True) # speeds up certain types of code
def simulate_trace(gamma = 1e-5,kappa = .01,fsim = 400,sim_points = 10240,seed = False):
    """Takes parameters and gives a simulated bead trajectory.
    .. math::
        dx = \frac{dt}{\alpha}\left(\sqrt{\frac{2\alpha k_BT}{dt}} F_L-\kappa x_{i-1}\right) \\
        x_i = x_{i-1}+dx
    """
    if seed != False:
        np.random.seed(seed)
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

def downsampled_trace(gamma = 1.e-5,kappa = .002,fsample = 100,N = 10240,seed=False):
    """
    Simulates a downsampled bead trajectory based on the overdamped Langevin equation.

    Parameters
    ----------
    gamma : float, optional
        Drag coefficient of the bead, by default 1.e-5
    kappa : float, optional
        Spring constant of the trap, by default .002
    fsample : int, optional
        Sampling frequency of the instrument, by default 100
    N : int, optional
        Number of points, by default 10240
    seed : int, optional
        Seed for random sampler, by default None

    Returns
    -------
    array
        Downsampled bead trajectory.
    """
    bin_size = 1000 # points to average when downsampling
    fsim = int(bin_size*fsample) # simulation frequency
    sim_points = N*bin_size # number of points to simulate
    xtrace = simulate_trace(gamma,kappa,fsim,sim_points,seed) # do simulation
    bins = xtrace.reshape((N,bin_size)) # reshape into bins 
    xtrace_ds = np.mean(bins,axis=1) # average bins to produce downsampled trace
    return xtrace_ds