import numpy as np
from numba import jit
@jit(nopython=True) # speeds up certain types of code
def simulate_trace(alpha = 1e-5,kappa = .002,fsim = 100,sim_points = 10240):
    """Takes parameters and gives a simulated probe trajectory."""
    dt = 1/fsim # time step between collecting probe position
    kT = 4.1    # thermal energy, in pN nm (approx. 4.1 at room temp.)
    F_L = np.random.standard_normal(sim_points) # Random Langevin force
    xtrace = np.zeros(sim_points) 
    for i in range(1,sim_points):
        dx =  dt/alpha * (np.sqrt(2*alpha*kT/dt)*F_L[i]-kappa*xtrace[i-1]) # Smoluchowski eqn
        xtrace[i] = xtrace[i-1] + dx
    return xtrace

def downsampled_trace(alpha = 1e-5,kappa = .002,fsample = 100,N = 10240):
    fsim = 3e5 # simulation frequency
    bin_size = np.floor(fsim/fsample).astype(int) # points to be averaged when downsampling
    sim_points = N*bin_size # number of points to simulate
    xtrace = simulate_trace(alpha,kappa,fsim,sim_points) # do simulation
    bins = xtrace.reshape((N,bin_size)) # reshape into bins 
    xtrace_ds = np.mean(bins,axis=1) # average bins to produce downsampled trace
    return xtrace_ds