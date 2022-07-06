import numpy as np

def ig(nu,tau,epsilon):
    """
    Calculates the inverse gamma
    """
    n = len(epsilon)
    nul = nu + n/2
    taul = tau + np.linalg.norm(epsilon)**2/2
    psi = 1/np.random.gamma(nul,1/taul)

    return psi