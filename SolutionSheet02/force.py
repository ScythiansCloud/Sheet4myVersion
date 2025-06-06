import settings
import numpy as np
from numba import njit, prange


@njit(parallel=True)
#@njit

#### unit of the force: (kcal/mole)/nm

def forceLJ(x, y, z, xlo, xhi, ylo, yhi, zlo, zhi, eps, sig, cutoff, epsWall, cutoffwall, sigWall):
    
    fx = np.zeros(shape=len(x))
    fy = np.zeros(shape=len(x))
    fz= np.zeros(shape=len(x))
    N = len(x)
    
    i = 0
    # sf2a = sig*sig / cutoff / cutoff
    # sf6a = sf2a * sf2a * sf2a
    #epotcut = 8.*settings.eps*sf6a*(sf6a - 1.)
    epot = 0
#    for i in range(N-1):
    for i in prange(N-1):
        j = i + 1
#        for j in range(i+1,N):
        for j in range(i+1, N):
            rijx = pbc(x[i], x[j],xlo, xhi)
            rijy = pbc(y[i], y[j],ylo, yhi)
            #rijz = pbc(z[i], z[j],zlo, zhi)
            
            rijz = z[j]-z[i]#pbc(z[i], z[j],zlo, zhi)
            
            r2 = rijx * rijx + rijy * rijy + rijz * rijz
            # calculate fx, fy, fz
            if r2 < cutoff * cutoff:
                sf2 = sig*sig / r2
                sf6 = sf2 * sf2 * sf2
                # epot += (8.*eps*sf6*(sf6 - 1.)) #-epotcut)
                ff = 48.*eps*sf6*(sf6 - 0.5)/r2
                fx[i] -= ff*rijx
                fy[i] -= ff*rijy
                fz[i] -= ff*rijz
                
                fx[j] += ff*rijx
                fy[j] += ff*rijy
                fz[j] += ff*rijz
    
    prefacW = 3*np.sqrt(3) *0.5 * epsWall
    sigfluidW3 = 3*sigWall**3
    sigfluidW9 = 9*sigWall**9

    fwall0 = 0 # initialize force on the wall (action = reaction)
    fwall2l = 0
    for i in prange(N):
        if z[i] < cutoffwall:
            f = prefacW* (sigfluidW9/z[i]**10 - sigfluidW3/z[i]**4)
            fz[i] += f
            fwall0 -= f
        if (zhi-cutoffwall)<z[i] < zhi:
            f = prefacW* (sigfluidW9/(zhi-z[i])**10 - sigfluidW3/(zhi-z[i])**4)
            fz[i] -= f
            fwall2l += f
        if zhi<z[i] or z[i]<0:
            print('particle escaped.....')
                
    return fx, fy, fz, fwall0, fwall2l, epot
  
@njit  
def pbc(xi, xj, xlo, xhi):
    
    l = xhi-xlo
    
    xi = xi % l
    xj = xj % l
    
    rij = xj - xi  
    if abs(rij) > 0.5*l:
        rij = rij - np.sign(rij) * l 
        
    return rij




# @njit
# def pbc(xi, xj, xlo, xhi):
#     l = xhi - xlo
#     rij = xi - xj
#     if rij > 0.5 * l:
#         rij -= l
#     elif rij < -0.5 * l:
#         rij += l
#     return rij
