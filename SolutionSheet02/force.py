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

    for i in prange(N):
        if z[i] < cutoffwall:
            fz[i] += prefacW* (sigfluidW9/z[i]**10 - sigfluidW3/z[i]**4)
        if (zhi-cutoff)<z[i] < zhi:
            fz[i] -= prefacW* (sigfluidW9/(zhi-z[i])**10 - sigfluidW3/(zhi-z[i])**4)
        if zhi<z[i] or z[i]<0:
            print('particle escaped.....')


    
    # for i in prange(N):
    #     d_bot  = z[i]        # z distance to left wall at z=0
    #     d_top = 2*zhi - z[i]  # z distance to right wall at z=2L

    #     if 0.0 < d_bot < cutoffwall:  # bottom
    #         # -dU/dd =  (9*sig9/d^10 - 3*sig3/d^4) * prefacW
    #         fz[i] += prefacW * (sigfluidW9/d_bot**10 - sigfluidW3/d_bot**4)
    #         epot  += prefacW * (sigfluidW9/d_bot**9  - sigfluidW3/d_bot**3)

    #     if 0.0 < d_top < cutoffwall:  # top
    #         # same dU/dd, but dd/dz = -1 so we subtract
    #         fz[i] -= prefacW * (sigfluidW9/d_top**10 - sigfluidW3/d_top**4)
    #         epot  += prefacW * (sigfluidW9/d_top**9  - sigfluidW3/d_top**3)
                
    return fx, fy, fz, epot
  
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
