import numpy as np
from ECRrectmodel import ECRrectmodel
from LikelihoodECRrect import LikelihoodECRrect
from Metropolis import METROPOLIS
from ig import ig
def MCMCECRig(ax, ay, az, z, t, kmin, kmax, k, SIGMAk, Dmin, Dmax, D, SIGMAD, N, nu, tau, thinfact):
    """

     The inputs of the function are:

     ax          = x dimension of the rectangle
     ay          = y dimension of the rectangle
     az          = z dimension of the rectangle
     Note: The code uses dimensionless parameters internally, so the user must
           use consistent units.
     z           = Data input as a function of t.
     t           = Time vector from data (independant value)
     kmin        = The smallest permissible k*
     kmax        = The largest permissible k*
     k           = Initiall value of k*
     SIGMAk      = Standard deviation of the draws on k
     Dmin        = The smallest permissible D*
     Dmax        = The largest permissible D*
     D           = Initiall value of D*
     SIGMAD      = Standard deviation of of the draws on D
     N           = Number of cycles to run the calibration (suggested = 100000.
                   There is a burn-in of around 5000, so the program returns
                   the results from the last N-5000) (This will take a lot of time)
     And it returns:

     kpost       = the posterior k* distribution
     Dpost       = the posterior D* distribution
     psipost     = the posterior psi distribution
     SUCCESSk    = The number of accepted k*
     SUCCESSD    = The number of accepted D*


    """
    print(' MCMCECRig is free software: you can redistribute it and/or modify')
    print(' it under the terms of the GNU General Public License as published by')
    print("the Free Software Foundation, either version 3 of the License, or")
    print(' (at your option) any later version.')

    #initialization and preallocation

    kpost = np.zeros((1, N))
    Dpost = np.zeros((1, N))
    psipost = np.zeros((1, N))
    likelies = np.zeros((1, N))

    kpost[0,0] = k
    Dpost[0,0] = D
    psipost[0,0] = tau/(nu+1)

    acceptedk = 0
    acceptedD = 0

    totacck = 0
    totaccD = 0

    y = ECRrectmodel(k, D, ax, ay, az, t)

    LOGLIKELY, EPSILON = LikelihoodECRrect(y, z, t, psipost[0][0])

    likelies[0,0] = LOGLIKELY

    for i in range(1,N):
        if np.mod(i,100) == 0:
            print('Currently on draw ' + str(i) + ' of ' + str(N) + ' total draws')
            print('Acceptace rate k = ' + str(acceptedk / 100))
            print('Acceptace rate D = ' + str(acceptedD / 100))

            if acceptedk/100 < 0.01:
                SIGMAk = SIGMAk/2
            elif acceptedk/100 > 0.1:
                SIGMAk = 1.5*SIGMAk
            if acceptedD/100 < 0.01:
                SIGMAD = SIGMAD/2
            elif acceptedD/100 > 0.1:
                SIGMAD = 1.5*SIGMAD
            totacck = totacck + acceptedk
            totaccD = totaccD + acceptedD
            acceptedD = 0
            acceptedk = 0
        proposedk = np.log10(kmin) - 1
        while proposedk < np.log10(kmin) or proposedk > np.log10(kmax):
            proposedk = np.random.normal(np.log10(kpost[0][i - 1]), SIGMAk)
        print(i)
        proposedk = 10**proposedk

        proposedy = ECRrectmodel(proposedk, Dpost[0][i - 1], ax, ay, az, t)

        proposedLikely, propsilon = LikelihoodECRrect(proposedy, z, t, psipost[0][i - 1])

        if METROPOLIS(proposedLikely, LOGLIKELY, thinfact):
            kpost[0][i] = proposedk
            LOGLIKELY = proposedLikely
            EPSILON = propsilon
            acceptedk = acceptedk + 1
        else:
            kpost[0][i] = kpost[0][i - 1]

        proposedD = np.log10(Dmin) - 1

        while proposedD < np.log10(Dmin) or proposedD > np.log10(Dmax):
            proposedD = np.random.normal(np.log10(Dpost[0][i - 1]), SIGMAD)

        proposedD = 10**proposedD

        proposedy = ECRrectmodel(kpost[0][i-1], proposedD, ax, ay, az, t)

        proposedLikely, propsilon = LikelihoodECRrect(proposedy, z, t, psipost[0][i - 1])

        if METROPOLIS(proposedLikely, LOGLIKELY, thinfact):
            Dpost[0][i] = proposedD
            LOGLIKELY = proposedLikely
            EPSILON = propsilon
            acceptedD = acceptedD + 1
        else:
            Dpost[0][i] = Dpost[0][i - 1]

        psipost[0][i] = ig(nu, tau, EPSILON)
        likelies[0][i] = LOGLIKELY


    successk = totacck/N
    successD = totaccD/N

    weights = (1 - thinfact)*likelies
    weights = weights - np.amax(weights)
    weights = np.exp(weights)

    return kpost, Dpost, psipost, weights, successk, successD