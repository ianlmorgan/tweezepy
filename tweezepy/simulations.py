import numpy as np
from numba import jit
@jit(nopython=True) # speeds up certain types of code
def simulate_trace(alpha = 1e-5,kappa = .001,fsim = 400,sim_points = 10240,seed = None):
    """Takes parameters and gives a simulated bead trajectory.
    .. math::
        dx = \frac{dt}{\alpha}\left(\sqrt{\frac{2\alpha k_BT}{dt}} F_L-\kappa x_{i-1}\right) \\
        x_i = x_{i-1}+dx
    """
    if fsim < kappa/alpha: # Simulation is unstable below this condition
        raise ValueError(r'The simulation frequency must be larger than kappa/alpha.')
    dt = 1/fsim # time step between collecting probe position
    kT = 4.1    # thermal energy, in pN nm (approx. 4.1 at room temp.)
    if seed is not None:
        np.random.seed(seed)
    F_L = np.random.standard_normal(sim_points-1) # Random Langevin force
    xtrace = np.zeros(sim_points) 
    for i in range(1,sim_points):
        dx =  dt/alpha * (np.sqrt(2*alpha*kT/dt)*F_L[i-1]-kappa*xtrace[i-1]) # Smoluchowski eqn
        xtrace[i] = xtrace[i-1] + dx
    return xtrace

def downsampled_trace(alpha = 1e-5,kappa = .002,fsample = 100,N = 10240,**kwargs):
    fsim = 3e5 # simulation frequency
    bin_size = np.floor(fsim/fsample).astype(int) # points to be averaged when downsampling
    sim_points = N*bin_size # number of points to simulate
    xtrace = simulate_trace(alpha,kappa,fsim,sim_points,**kwargs) # do simulation
    bins = xtrace.reshape((N,bin_size)) # reshape into bins 
    xtrace_ds = np.mean(bins,axis=1) # average bins to produce downsampled trace
    return xtrace_ds