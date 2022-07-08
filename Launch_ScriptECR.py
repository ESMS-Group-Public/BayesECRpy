import numpy as np
from ECRrectmodel import ECRrectmodel
from MCMCECRig import MCMCECRig
import matplotlib.pyplot as py

Matrix = np.loadtxt('Graphically_Extracted_dataREF_24.csv',dtype=float,delimiter=',')
t3 = Matrix[:,0]
z = Matrix[:,1]

#Prior
ax= 4
ay= .15
az= .15


kmin     =   1e-9
kmax     =   1
k        =   1e-4
SIGMAk   =   0.1

Dmin     =   1e-10
Dmax     =   1e-2
D        =   1e-5
SIGMAD   =   0.1

ps      = 0.03
N        = 10000
nu       = 1000
tau      = 1000*ps**2

thinfactor = .05

kpost, Dpost, psipost, weights, successk, successD = MCMCECRig(ax, ay, az, z, t3, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfactor)

covdraws = 50
drawind=np.random.randint(low=1, high=N-1000,size=covdraws)

py.plot(t3,z, 'r')

for i in range(1,covdraws):
    y = ECRrectmodel(kpost[0][1000 + drawind[i]], Dpost[0][1000 + drawind[i]], ax, ay, az, t3)
    py.plot(t3, y)

py.show()

#Scatter plot

scatdraws = 1000
scatind = np.random.randint(low=1, high=N-1000,size=scatdraws)
scattervecs = np.zeros((scatdraws,2))
for i in range(1, scatdraws):
    scattervecs[i, :] = [kpost[0][scatind[i] + 1000], Dpost[0][scatind[i] + 1000]]

py.scatter(scattervecs[:, 0], scattervecs[:, 1],alpha=.2)

py.show()