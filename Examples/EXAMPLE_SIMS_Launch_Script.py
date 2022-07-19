import numpy as np
from MCMCSIMSig import MCMCSIMSig
# Example of SIMS use cited as REF 23
# Original paper htat details experiment can be found at
# https://onlinelibrary.wiley.com/doi/abs/10.1002/fuce.201300087

# A dataset is included for this tutorial. It is pH2O-dependent 1.txt. This SIMS data set was produced by Dr. Roger
# De-Souze that has already been analyzed
# Tutorial is based on pH2O = 220mbar experiment

#Importing Data as row vectors
H220x = np.loadtxt('H220x.txt', dtype=float, delimiter=',')
H220z = np.loadtxt('H220z.txt', dtype=float, delimiter=',')

# First row is set of sample thicknesses in cm
# Second row is the set of O18 fractions

H220t = 5640

# Now all that is needed is to enter the values we want to use for the other
# inputs.  These can be entered either in the workspace or in the command
# window.

# Note the wide bounds on k and D
# k, surface reaction rate
kmin     = 1e-15      # Smallest permissible k*
kmax     = 1e-2       # Largest permissible k*
k        = 1e-6       # Initial value of k*
SIGMAk   = .1         # Standard deviation of the draws on k
# D, bulk diffusion constant
Dmin     = 1e-15      # Smallest permissible D*
Dmax     = 1e-2       # Largest permissible D*
D        = 5e-6       # Initial value of D*
SIGMAD   = .1         # Standard deviation of the draws on D

ps       = 0.03       # Standard Deviation of Observational error
N        = 5000       # Number of cycles to run the calibration
nu       = 1000       # Shape parameter of the inverse gamma distribution
tau      = 1001*ps**2 # Scale parameter of the inverse gamma distribution
thinfact = 1          # Thinning factor b/w [0,1]: fraction of data to be used for calibration

MCMCSIMSig(H220x, H220z, H220t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact)