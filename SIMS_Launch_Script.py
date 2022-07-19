import numpy as np
from MCMCSIMSig import MCMCSIMSig

# LAUNCHER FOR SIMS
# Please see Readme and Example case of SIMS for more information


#Import data as row vectors [1xN]
x        =     # Depth in cm
z        =     # tracer site fraction, distance corresponds to x

t        =     # time associated with the SIMS dataset [1x1]

# Priors, or initial guesses at values

# k, surface reaction rate
kmin     =
kmax     =
k        =
SIGMAk   =
# D, bulk diffusion constant
Dmin     =
Dmax     =
D        =
SIGMAD   =

# SAMPLE OF VALUES FOR IG DISTRIBUTION; Experiment with your own values
ps       =            # StD of Obs. Variance
N        = 5000       # Number of cycles to run the calibration
nu       = 4          # Shape parameter of the ig distribution
tau      = 5*ps**2    # Scale parameter of the ig distribution

thinfact =            # Thinning factor between [0,1]

MCMCSIMSig(H220x, H220z, H220t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact)