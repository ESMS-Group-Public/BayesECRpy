import numpy as np

def conf(k,D,L,psi,P):

    burnin_drop = 4999
    k[np.arange(1,burnin_drop,1)] = []
    D[np.arange(1,burnin_drop,1)] = []
    L[np.arange(1,burnin_drop,1)] = []
    psi[np.arange(1,burnin_drop,1)] = []
    Lsort = np.sort(L)
    a = np.length(L)
    b = np.floor((100 - P) / 100 * a)
    c = Lsort[b]
    j = 1
    k95 = 0
    D95 = 0
    psi95 = 0
    D95t = 0
    for i in range(1,np.length(L)):
        if L[i] > c:
            k95[i] = k[i]
            D95[i] = D[i]
            psi95[i] = psi[i]
        elif k95[i] == 0:
            D95[i] = 0
            psi[i] = 0
    k = 1
    p = 1
    while k <= np.length(D95):
        if D95[k] != 0:
            D95t[p] = D95[k]
            p = p + 1
        k = k+1

    k = 1
    p = 1
    psi95t = 0
    while k <= np.length(psi95):
        if psi95[k] != 0:
            psi95t[p] = psi95[k]
            p = p + 1
        k = k+1


    k = 1
    p = 1
    k95t = 0
    while k <= np.length(k95):
        if k95[k] != 0:
            k95t[p] = k95[k]
            p = p + 1
        k = k+1

    D95 = D95t
    k95 = k95t
    psi95 = psi95t

    return k95, D95, psi95