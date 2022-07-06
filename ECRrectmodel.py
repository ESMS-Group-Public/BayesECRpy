import numpy as np
import scipy.optimize as op
def ECRrectmodel(k, D, ax, ay, az, t):
    """
    Solves the ECR rectangle model as described in Yasuda and Hishinuma's
    1996 paper in J. Solid State Chem
    Returns a vector y = M(t)/M(inf)
    k = effective exchange coefficient
    D = Effective Diffusion coefficient
    ax, ay, az are the x, y, and z dimensions of the rectangle
    t = is the time vector
    """
    Lx = ax*k/D

    TERMx = 1
    SUMx = 1

    def funcx(bn):
        return bn*np.tan(bn)-Lx

    def funcy(bn):
        return bn*np.tan(bn)-Ly

    def funcz(bn):
        return bn*np.tan(bn)-Lz

    box = op.root_scalar(funcx,bracket = [-.00000000001, (np.pi/2)-.000000001], method='toms748') #fzero func
    box = box.root
    TERMx = 2*Lx**2*np.exp(-(box**2)*D*(t)/(ax**2))/(box**2*(box**2+Lx**2+Lx))
    SUMx = TERMx

    ix = 2


    while np.amax(TERMx/SUMx) > .0001:
        bn = op.root_scalar(funcx,bracket = [(ix-1)*np.pi-.00000000001, (2*(ix-1)+1)*(np.pi/2)-.000000001], method='toms748') #fzero func #FSolve
        bn = bn.root
        TERMx = 2*Lx**2*np.exp(-(bn**2)*D*(t)/(ax**2))/(bn**2*(bn**2+Lx**2+Lx))
        SUMx = SUMx + TERMx
        ix = ix + 1

    Ly = ay*k/D

    TERMy = 1
    SUMy = 1

    boy = op.root_scalar(funcy,bracket = [-.00000000001, (np.pi/2)-.000000001], method='toms748') #fzero func
    boy = boy.root
    TERMy  = 2*Ly**2*np.exp(-(boy**2)*D*(t)/(ay**2))/(boy**2*(boy**2+Ly**2+Ly))

    SUMy = TERMy

    iy = 2
    while np.amax(TERMy/SUMy) > .0001:
        bny = op.root_scalar(funcy,bracket = [(iy-1)*np.pi-.00000000001, (2*(iy-1)+1)*(np.pi/2)-.000000001], method='toms748') #fzero func
        bny = bny.root
        TERMy = 2*Ly**2*np.exp(-(bny**2)*D*(t)/(ay**2))/(bny**2*(bny**2+Ly**2+Ly))
        SUMy = SUMy + TERMy
        iy = iy + 1

    Lz = az*k/D

    TERMz = 1
    SUMz = 1

    boz = op.root_scalar(funcz,bracket = [-.00000000001, (np.pi/2)-.000000001], method='toms748') #fzero func
    boz = boz.root
    TERMz = 2*Lz**2*np.exp(-(boz**2)*D*(t)/(az**2))/(boz**2*(boz**2+Lz**2+Lz))

    SUMz = TERMz

    iz = 2

    while np.amax(TERMz/SUMz) > .0001:
        bnz = op.root_scalar(funcz,bracket = [(iz-1)*np.pi-.00000000001, (2*(iz-1)+1)*(np.pi/2)-.000000001], method='toms748') #fzero func
        bnz = bnz.root
        TERMz = 2*Lz**2*np.exp(-(bnz**2)*D*(t)/(az**2))/(bnz**2*(bnz**2+Lz**2+Lz)) #matrix comp
        SUMz = SUMz + TERMz
        iz = iz + 1;

    SUM3D = SUMx*SUMy*SUMz #matrix comp
    y = 1-SUM3D

    return y


