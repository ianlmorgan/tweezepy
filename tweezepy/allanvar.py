"""
This module performs operations associated with calculating the Allan variance. 

"""

import numpy as np
import scipy
import warnings 

def m_generator(N,taus = 'octave'):
    assert taus in ['all','octave','decade'], "taus must be either all, octave, or decade."
    if taus == 'all':
        # all-tau sampling not particularly useful but why not?
        maxn = N//2
        m = np.linspace(1.0,maxn,maxn,dtype='int')
    elif taus == 'octave':
        # octave sampling break bin sizes, m, into powers of 2^n
        maxn = int(np.floor(np.log2(N/2))) # m =< N/2
        m = np.logspace(0,maxn-1,maxn,base=2,dtype='int')  #bin sizes
    elif taus == 'decade':
        # again not particularly useful, but why not?
        maxn = int(np.floor(np.log10(N/2)))
        m = [np.array([1,2,4])*k for k in np.logspace(0,maxn,maxn+1,base=10,dtype='int')]
        m = np.ravel(m)
    return m

def calc_avar(phase,rate,mj,step):
    d2 = phase[2 * mj::step]
    d1 = phase[1 * mj::step]
    d0 = phase[::step]
    nmin = min(len(d0), len(d1), len(d2))
    v_arr = d2[:nmin] - 2 * d1[:nmin] + d0[:nmin]
    avar = np.sum(v_arr*v_arr)/(2.*nmin*pow(mj/rate,2))
    return avar

def avar(data,rate = 1.0,taus = 'octave', overlapping = True, edf = 'approx'):
    """
    Calculate standard allan variance 
    Takes an array of bead positions.
    Returns the taus, edfs, and oavs.

    Parameters
    ----------
    data : array, series, or list
        1-D array of numbers.
    rate : float
        Frequency of acquisition.

    Returns
    -------
    taus : array
        Observation times.
    edfs : array
        Equivalent degrees of freedom.
    oavs : array
        Allan variance.
    
    Notes
    -----
    Adapted from Allantools
    """
    assert type(overlapping) == bool, 'overlapping keyword argument should be a boolean.'
    assert edf in ['approx','real'], 'edf keyword argument should be approx or real.'
    rate = float(rate)
    data = np.asarray(data) # convert to numpy array
    N = len(data)
    m = m_generator(N,taus = taus)
    n = N - 2*m + 1 
    m = m[n>=2]
    taus = m/rate # tau = m*tau_c
    if edf == 'real':
        edfs = np.empty(len(m))
        for i,mj in enumerate(m):
            if N//mj > 32:
                alpha_int = noise_id(data,mj)[0]
            if (alpha_int<3) and (alpha_int>-3):
                edfs[i] = edf_greenhall(alpha_int,2,mj,N)
            else:
                warnings.warn('Real edf failed to identify noise for %s. Falling back to approximate edf.'%mj)
                edfs[i] = edf_approx(N,mj)
    elif edf == 'approx':
        edfs = edf_approx(N,m)
    else:
        print('edf keyword argument %s not recognized.'%edf)
        raise UserWarning
    # Calculate phasedata from Eq. 18b (in erratum)
    phase = np.cumsum(data)/rate  # integrate positions, converting frequency to phase data
    phase = np.insert(phase, 0, 0) # phase data should start at 0
    assert len(phase) > 0, "Data array length is too small: %i" % len(phase)
    # Calculate oav from Eq. 18a (in erratum)
    if overlapping:
        oavs = np.array([calc_avar(phase,rate,mj,1) for mj in m])
    elif not overlapping:
        oavs = np.array([calc_avar(phase,rate,mj,mj) for mj in m])
    return taus,edfs,oavs

def totvar(data, rate=1.0, taus='octave',edf = 'approx'):
    """ Total variance.
        Better confidence at long averages for Allan variance.
    .. math::
        \\sigma^2_{TOTDEV}( m\\tau_0 ) = { 1 \\over 2 (m\\tau_0)^2 (N-2) }
            \\sum_{i=2}^{N-1} ( {x}^*_{i-m} - 2x^*_{i} + x^*_{i+m} )^2
    where :math:`x^*_i` is a new time-series of length :math:`3N-4`
    derived from the original phase time-series :math:`x_n` of
    length :math:`N` by reflection at both ends.
    Parameters
    ----------
    data: np.array
        Bead positions.
    rate: float
        The sampling rate for phase or frequency, in Hz
    taus: np.array
        Array of tau values for which to compute measurement
    """
    rate = float(rate)
    data = np.asarray(data) # make sure data is an array, not a series or list
    m = m_generator(len(data),taus=taus)
    taus = m/rate

    phase = np.cumsum(data)/rate  # convert frequency data to phase data
    phase = np.insert(phase, 0, 0) # phase data should start at 0
    N = len(phase)
    assert N > 0, "Data array length is too small: %i" % len(phase)
    
    # totvar requires a new dataset
    # Begin by adding reflected data before dataset
    x1 = 2.0 * phase[0] * np.ones((N - 2,))
    x1 = x1 - phase[1:-1]
    x1 = x1[::-1]

    # Reflected data at end of dataset
    x2 = 2.0 * phase[-1] * np.ones((N - 2,))
    x2 = x2 - phase[1:-1][::-1]

    # check length of new dataset
    assert len(x1)+len(phase)+len(x2) == 3*N - 4

    # Combine into a single array
    x = np.zeros((3*N - 4))
    x[0:N-2] = x1
    x[N-2:2*(N-2)+2] = phase  # original data in the middle
    x[2*(N-2)+2:] = x2
    
    mid = len(x1)
    tvars= np.array([calc_totvar(x,rate,N,mid,mj) for mj in m])
    if edf == 'approx':
        edfs = edf_approx(N,m)
    elif edf == 'real':
        edfs = np.empty(len(m))
        for idx,mj in enumerate(m):
            if N//mj >= 32:
                alpha_int = noise_id(data,mj)[0]
            tvars[idx] /= 1-totvar_bias(alpha_int)*(mj/N)
            if (alpha_int <= 0) and (alpha_int >= -2):
                edfs[idx] = edf_totdev(N,mj,alpha_int)
            else:
                warnings.warn('Real edf failed for %s.'%mj)
                edfs[idx] = edf_approx(N,mj)
    return taus,edfs,tvars
def totvar_bias(alpha_int):
    if alpha_int == -2:
        return .75
    elif alpha_int == -1:
        return  0.481
    else:
        return 0
def calc_totvar(x,rate,N,mid,mj):
    d0 = x[mid + 1:]
    d1 = x[mid + mj + 1:]
    d1n = x[mid - mj + 1:]
    nmin = min(len(d0), len(d1), len(d1n))
    v_arr = d1n[:nmin] - 2.0 * d0[:nmin] + d1[:nmin]
    var = np.sum(v_arr[:mid] * v_arr[:mid])
    var /= float(2 * pow(mj / rate, 2) * (N-2))
    return var

########################################################################
# Noise Identification Lag1 autorrelation

def noise_id(x,af, dmin = 0, dmax = 2):
    """
    Lag1 autocorrelation algorithm for determining powerlaw noise type.

    Parameters
    ----------
    x : numpy.array
        trace
    af : int
        averaging factor
    dmin : int, optional
        minimum number of differentiations, by default 0
    dmax : int, optional
        maximum number of differentiations, by default 2

    Returns
    -------
    (alpha_int,alpha) : tuple
        [description]
    alpha_int : int
        integer power-law noise
    alpha : float
        float power-law noise

    """
    # Split time series into average positions of nonoverlapping bins
    N = len(x)
    x = x[:N//af * af].reshape((N//af,af)) # cut to correct lengths
    x = np.average(x,axis=1)
    # require minimum length for time-series
    if len(x) < 32:
        #print(("noise_id() Can't determine noise-ID for"
        #       " time-series length= %d") % len(x))
        return np.nan,np.nan
        #raise NotImplementedError
    d = 0 # number of differentiations
    while True:
        r1 = np.corrcoef(x[1:],x[:-1])[0][1]
        rho = r1/(1.0+r1)
        if d >= dmin and (rho < 0.25 or d >= dmax):
            alpha = -2.*(rho+d)
            alpha_int = int(-1.0*np.round(2*rho) - 2.0*d)      
            return alpha_int,alpha
            
        else:
            x = np.diff(x)
            d = d + 1
            
def edf_greenhall(alpha, d, m, N,
                  overlapping=False, modified=False, verbose=False):
    """ returns Equivalent degrees of freedom - couresy of allantools
        Parameters
        ----------
        alpha: int
            noise type, +2...-4
        d: int
            1 first-difference variance
            2 Allan variance
            3 Hadamard variance
            require alpha+2*d>1
        m: int
            averaging factor
            tau = m*tau0 = m*(1/rate)
        N: int
            number of phase observations (length of time-series)
        overlapping: bool
            True for oadev, ohdev
        modified: bool
            True for mdev, tdev
        Returns
        -------
        edf: float
            Equivalent degrees of freedom
        Greenhall, Riley, 2004
        https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20050061319.pdf
        UNCERTAINTY OF STABILITY VARIANCES BASED ON FINITE DIFFERENCES
        Notes
        -----
        Based on allantools https://github.com/aewallin/allantools

        Used for the following deviations
        (see http://www.wriley.com/CI2.pdf page 8)
        adev()
        oadev()
        mdev()
        tdev()
        hdev()
        ohdev()
    """

    if modified:
        F = 1  # F filter factor, 1 modified variance, m unmodified variance
    else:
        F = int(m)
    if overlapping:  # S stride factor, 1 nonoverlapped estimator,
        S = int(m)  # m overlapped estimator (estimator stride = tau/S )
    else:
        S = 1
    assert(alpha+2*d > 1.0)
    L = m/F+m*d  # length of filter applied to phase samples
    M = 1 + np.floor(S*(N-L) / m)
    J = min(M, (d+1)*S)
    J_max = 100
    r = M/S
    if int(F) == 1 and modified:  # case 1, modified variances, all alpha
        if J <= J_max:
            inv_edf = (1.0/(pow(greenhall_sz(0, 1, alpha, d), 2)*M)) * \
                       greenhall_BasicSum(J, M, S, 1, alpha, d)
            if verbose:
                print("case 1.1 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        elif r > d+1:
            (a0, a1) = greenhall_table1(alpha, d)
            inv_edf = (1.0/r)*(a0-a1/r)
            if verbose:
                print("case 1.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            inv_edf = ((1.0/(pow(greenhall_sz(0, F, alpha, d), 2)*J_max)) *
                       greenhall_BasicSum(J_max, J_max, m_prime, 1, alpha, d))
            if verbose:
                print("case 1.3 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
    elif int(F) == int(m) and int(alpha) <= 0 and not modified:
        # case 2, unmodified variances, alpha <= 0
        if J <= J_max:
            if m*(d+1) <= J_max:
                m_prime = m
                variant = "a"
            else:
                m_prime = float('inf')
                variant = "b"

            inv_edf = ((1.0/(pow(greenhall_sz(0, m_prime, alpha, d), 2)*M)) *
                       greenhall_BasicSum(J, M, S, m_prime, alpha, d))
            if verbose:
                print("case 2.1%s edf= %3f" % (variant, float(1.0/inv_edf)))
            return 1.0/inv_edf
        elif r > d+1:
            (a0, a1) = greenhall_table2(alpha, d)
            inv_edf = (1.0/r)*(a0-a1/r)
            if verbose:
                print("case 2.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            inv_edf = (
                (1.0/(pow(greenhall_sz(0, float('inf'), alpha, d), 2)*J_max)) *
                greenhall_BasicSum(
                    J_max, J_max, m_prime, float('inf'), alpha, d))
            if verbose:
                print("case 2.3 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
    elif int(F) == int(m) and int(alpha) == 1 and not modified:
        # case 3, unmodified variances, alpha=1
        if J <= J_max:
            # note: m<1e6 to avoid roundoff
            inv_edf = ((1.0/(pow(greenhall_sz(0, m, 1, d), 2)*M)) *
                       greenhall_BasicSum(J, M, S, m, 1, d))
            if verbose:
                print("case 3.1 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        elif r > d+1:
            (a0, a1) = greenhall_table2(alpha, d)
            (b0, b1) = greenhall_table3(alpha, d)
            inv_edf = (1.0/(pow(b0+b1*np.log(m), 2)*r))*(a0-a1/r)
            if verbose:
                print("case 3.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        else:
            m_prime = J_max/r
            (b0, b1) = greenhall_table3(alpha, d)
            inv_edf = (
                (1.0/(pow(b0+b1*np.log(m), 2)*J_max)) *
                greenhall_BasicSum(J_max, J_max, m_prime, m_prime, 1, d))
            if verbose:
                print("case 3.3 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
    elif int(F) == int(m) and int(alpha) == 2 and not modified:
        # case 4, unmodified variances, alpha=2
        K = np.ceil(r).astype('int')
        if K <= d:
            inv_edf = (1.0 + 
                      (2.0/scipy.special.binom(2*d,d)) * 
                      np.sum([(1.0-k/r)*pow(scipy.special.binom(2*d,d-k),2) for k in range(1,K)]))
            if verbose:
                print("case 4.1 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf
        else:
            a0 = (scipy.special.binom(4*d, 2*d) /
                  pow(scipy.special.binom(2*d, d), 2))
            a1 = d/2.0
            inv_edf = (1.0/M)*(a0-a1/r)
            if verbose:
                print("case 4.2 edf= %3f" % float(1.0/inv_edf))
            return 1.0/inv_edf

    print("greenhall_edf() no matching case!")
    raise NotImplementedError

def greenhall_BasicSum(J, M, S, F, alpha, d):
    """ Eqn (10) from Greenhall2004 """
    first = pow(greenhall_sz(0, F, alpha, d), 2)
    second = ((1-float(J)/float(M)) *
              pow(greenhall_sz(float(J)/float(S), F, alpha, d), 2))
    third = 0
    for j in range(1, int(J)):
        third += (2 * (1.0-float(j)/float(M)) *
                  pow(greenhall_sz(float(j) / float(S), F, alpha, d), 2))
    return first+second+third

def greenhall_sz(t, F, alpha, d):
    """ Eqn (9) from Greenhall2004 """
    if d == 1:
        a = 2*greenhall_sx(t, F, alpha)
        b = greenhall_sx(t-1.0, F, alpha)
        c = greenhall_sx(t+1.0, F, alpha)
        return a-b-c
    elif d == 2:
        a = 6*greenhall_sx(t, F, alpha)
        b = 4*greenhall_sx(t-1.0, F, alpha)
        c = 4*greenhall_sx(t+1.0, F, alpha)
        dd = greenhall_sx(t-2.0, F, alpha)
        e = greenhall_sx(t+2.0, F, alpha)
        return a-b-c+dd+e
    elif d == 3:
        a = 20.0*greenhall_sx(t, F, alpha)
        b = 15.0*greenhall_sx(t-1.0, F, alpha)
        c = 15.0*greenhall_sx(t+1.0, F, alpha)
        dd = 6.0*greenhall_sx(t-2.0, F, alpha)
        e = 6.0*greenhall_sx(t+2.0, F, alpha)
        f = greenhall_sx(t-3.0, F, alpha)
        g = greenhall_sx(t+3.0, F, alpha)
        return a-b-c+dd+e-f-g

    assert(0)  # ERROR

def greenhall_sx(t, F, alpha):
    """ Eqn (8) from Greenhall2004
    """
    if F == float('inf'):
        return greenhall_sw(t, alpha+2)
    a = 2*greenhall_sw(t, alpha)
    b = greenhall_sw(t-1.0/float(F), alpha)
    c = greenhall_sw(t+1.0/float(F), alpha)

    return pow(F, 2)*(a-b-c)


def greenhall_sw(t, alpha):
    """ Eqn (7) from Greenhall2004
    """
    alpha = int(alpha)
    if alpha == 2:
        return -np.abs(t)
    elif alpha == 1:
        if t == 0:
            return 0
        else:
            return pow(t, 2)*np.log(np.abs(t))
    elif alpha == 0:
        return np.abs(pow(t, 3))
    elif alpha == -1:
        if t == 0:
            return 0
        else:
            return pow(t, 4)*np.log(np.abs(t))
    elif alpha == -2:
        return np.abs(pow(t, 5))
    elif alpha == -3:
        if t == 0:
            return 0
        else:
            return pow(t, 6)*np.log(np.abs(t))
    elif alpha == -4:
        return np.abs(pow(t, 7))

    assert(0)  # ERROR


def greenhall_table3(alpha, d):
    """ Table 3 from Greenhall 2004 """
    assert(alpha == 1)
    idx = d-1
    table3 = [(6.0, 4.0), (15.23, 12.0), (47.8, 40.0)]
    return table3[idx]


def greenhall_table2(alpha, d):
    """ Table 2 from Greenhall 2004 """
    row_idx = int(-alpha+2)  # map 2-> row0 and -4-> row6
    assert(row_idx in [0, 1, 2, 3, 4, 5])
    col_idx = int(d-1)
    table2 = [
        # alpha = +2:
        [(3.0/2.0, 1.0/2.0), (35.0/18.0, 1.0), (231.0/100.0, 3.0/2.0)],
        # alpha = +1:
        [(78.6, 25.2), (790.0, 410.0), (9950.0, 6520.0)],
        # alpha = 0:
        [(2.0/3.0, 1.0/6.0), (2.0/3.0, 1.0/3.0), (7.0/9.0, 1.0/2.0)],
        # alpha = -1:
        [(-1, -1), (0.852, 0.375), (0.997, 0.617)],
        # alpha = -2
        [(-1, -1), (1.079, 0.368), (1.033, 0.607)],
        # alpha = -3
        [(-1, -1), (-1, -1), (1.053, 0.553)],
        # alpha = -4
        [(-1, -1), (-1, -1), (1.302, 0.535)]]
    # print("table2 = ", table2[row_idx][col_idx])
    return table2[row_idx][col_idx]

def greenhall_table1(alpha, d):
    """ Table 1 from Greenhall 2004 """
    row_idx = int(-alpha+2)  # map 2-> row0 and -4-> row6
    col_idx = int(d-1)
    table1 = [
        # alpha = +2:
        [(2.0/3.0, 1.0/3.0), (7.0/9.0, 1.0/2.0), (22.0/25.0, 2.0/3.0)],
        # alpha = +1:
        [(0.840, 0.345), (0.997, 0.616), (1.141, 0.843)],
        # alpha = 0:
        [(1.079, 0.368), (1.033, 0.607), (1.184, 0.848)],
        # alpha = -1:
        [(-1, -1), (1.048, 0.534), (1.180, 0.816)],
        # alpha = -2
        [(-1, -1), (1.302, 0.535), (1.175, 0.777)],
        # alpha = -3
        [(-1, -1), (-1, -1), (1.194, 0.703)],
        # alpha = -4
        [(-1, -1), (-1, -1), (1.489, 0.702)]]
    # print("table1 = ", table1[row_idx][col_idx])
    return table1[row_idx][col_idx]

def edf_totdev(N, m, alpha):
    """ Equivalent degrees of freedom for Total Deviation
        NIST SP1065 page 41, Table 7
    """
    alpha = int(alpha)
    if alpha in [0, -1, -2]:
        # alpha  0 WFM
        # alpha -1 FFM
        # alpha -2 RWFM
        NIST_SP1065_table7 = [(1.50, 0.0), (1.17, 0.22), (0.93, 0.36)]
        (b, c) = NIST_SP1065_table7[int(abs(alpha))]
        return b*(float(N)/float(m))-c
    # alpha outside 0, -1, -2:
    return edf_simple(N, m, alpha)

def edf_simple(N, m, alpha, pedantic = False):
    """Equivalent degrees of freedom.
    Simple approximate formulae.
    Parameters
    ----------
    N : int
        the number of phase samples
    m : int
        averaging factor, tau = m * tau0
    alpha: int
        exponent of f for the frequency PSD:
        'wp' returns white phase noise.             alpha=+2
        'wf' returns white frequency noise.         alpha= 0
        'fp' returns flicker phase noise.           alpha=+1
        'ff' returns flicker frequency noise.       alpha=-1
        'rf' returns random walk frequency noise.   alpha=-2
        If the input is not recognized, it defaults to idealized, uncorrelated
        noise with (N-1) degrees of freedom.
    Notes
    -----
       S. Stein, Frequency and Time - Their Measurement and
       Characterization. Precision Frequency Control Vol 2, 1985, pp 191-416.
       http://tf.boulder.nist.gov/general/pdf/666.pdf
       
       Modified from allantools.
    Returns
    -------
    edf : float
        Equivalent degrees of freedom
    """

    N = float(N)
    m = float(m)
    if alpha in [2, 1, 0, -1, -2]:
        # NIST SP 1065, Table 5
        if alpha == +2:
            edf = (N + 1) * (N - 2*m) / (2 * (N - m))

        if alpha == 0:
            edf = (((3 * (N - 1) / (2 * m)) - (2 * (N - 2) / N)) *
                   ((4*pow(m, 2)) / ((4*pow(m, 2)) + 5)))

        if alpha == 1:
            a = (N - 1)/(2 * m)
            b = (2 * m + 1) * (N - 1) / 4
            edf = np.exp(np.sqrt(np.log(a) * np.log(b)))

        if alpha == -1:
            if m == 1:
                edf = 2 * (N - 2)/(2.3 * N - 4.9)
            if m >= 2:
                edf = 5 * N**2 / (4 * m * (N + (3 * m)))

        if alpha == -2:
            a = (N - 2) / (m * (N - 3)**2)
            b = (N - 1)**2
            c = 3 * m * (N - 1)
            d = 4 * m**2
            edf = a * (b - c + d)

    else:
        edf = (N/m - 1) # assume correlated noise
        if pedantic == True:
            print("Noise type not recognized."
                  " Defaulting to N/m - 1 degrees of freedom.")

    return edf


def edf_approx(N,mj):
    return N//mj-1
