import numpy as np

def LikelihoodECRrect(y, z, t, psi):
    """
        This function calculates the likelihood given the
        results from ECRrectmodel
    """
    p = len(t)

    # Check size
    epsilon = z - y
    LOGLIKELY = -p/2*np.log(psi)-(1/(2*psi))*np.dot(epsilon.reshape(len(epsilon),1).T,epsilon.reshape(len(epsilon),1))

    return LOGLIKELY, epsilon