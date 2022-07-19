import numpy as np
from ECRrectmodel import ECRrectmodel
from MCMCECRig import MCMCECRig
import matplotlib.pyplot as py

# Example of ECR use cited as REF 25
# Original paper that details experiment can be found at
# https://www.sciencedirect.com/science/article/abs/pii/S0167273800005543

# Data Extraction
# Data is formatted to two row vectors
Matrix = np.loadtxt('Graphically_Extracted_data[REF_24].csv',dtype=float,delimiter=',')
t = np.array(Matrix[:,0]).T
z = np.array(Matrix[:,1]).T

# Dimensions of rectangular specimen, in this case LSCF
# Use consistent units
ax = 4
ay = .15
az = .15

# Priors, or initial guesses at values
# since minimum values are small, guesses can be an order of magnitude
# k, surface reaction rate
kmin     =   1e-9       # Smallest permissible k*
kmax     =   1          # Largest permissible k*
k        =   1e-4       # Initial value of k*
SIGMAk   =   0.1        # Standard Deviation of the draws on k
# D, bulk diffusion constant
Dmin     =   1e-10      # Smallest permissible D*
Dmax     =   1e-2       # Largest permissible D*
D        =   1e-5       # Initial value of D*
SIGMAD   =   0.1        # Standard Deviation of the draws on D

ps      = 0.02          # Standard Deviation of Observational error
N        = 5000         # Number of cycles to run the calibration
nu       = 1000         # Shape parameter of the inverse gamma prior
tau      = 1010*ps**2   # Scale parameter of the inverse gamma prior

thinfact = .15          # Thinning factor b/w [0,1]: fraction of data to be used for calibration

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