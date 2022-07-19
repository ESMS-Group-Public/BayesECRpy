import numpy as np
from ECRrectmodel import ECRrectmodel
from MCMCECRig import MCMCECRig
import matplotlib.pyplot as py

# LAUNCHER FOR ECR
# Please see Readme and Example case of ECR for more informtaion

#Reduce your data to two row vectors
t        =      # Time [1xN]
z        =      # Row vector as a function of t [1xN]

# Dimensions of rectangular specimen
# Use consistent units
ax       =
ay       =
az       =

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

#SAMPLE OF VALUES FOR UG DISTRIBUTION; Experiment with your own values
ps       = 0.02         # StD of Obs. Variance
N        = 5000         # Number of cycles to run the calibration
nu       = 1000         # Shape parameter of the ig prior
tau      = 1010*ps**2   # Scale parameter of the ig prior

thinfact =              # Thinning factor b/w [0,1]

kpost, Dpost, psipost, weights, successk, successD = MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact)

# Plots
covdraws = 50
drawind=np.random.randint(low=1, high=N-1000,size=50)

py.plot(t,z, 'r.')

for i in range(1,covdraws):
    y = ECRrectmodel(kpost[0][1000 + drawind[i]], Dpost[0][1000 + drawind[i]], ax, ay, az, t)
    py.plot(t, y)

py.show()

# Scatter plot

scatdraws = 1000
scatind = np.random.randint(low=1, high=N-1000,size=scatdraws)
scattervecs = np.zeros((scatdraws,2))
for i in range(1, scatdraws):
    scattervecs[i, :] = [kpost[0][scatind[i] + 1000], Dpost[0][scatind[i] + 1000]]

py.scatter(scattervecs[:, 0], scattervecs[:, 1],alpha=.2)

py.show()