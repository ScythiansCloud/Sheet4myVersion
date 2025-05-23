# module containing integration sheme:
"""
     --------- units --------- 
     [r] = nm; [t] = fs; [epsilon] = kcal/mole; [m] = gram/mole;  [F] = (kcal/mole)/nm; 
     [T] = K; [v] = nm/fs; [k_B] = (kcal/mole)/K
     conversion factor:
         from kcal*fs*fs/gram/nm to nm: 4.1868e-06
         from kcal*fs/gram/nm to nm/fs: 4.1868e-06
"""

import settings
import force
import numpy as np
from numba import njit, prange

@njit(parallel=True)
def VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, xlo, xhi, ylo, yhi, zlo, zhi, eps, sigma, cutoff, deltat, mass, epsWall, cutoffwall, sigWall):

    # conversion factor
    fx0 = np.zeros(shape=len(x))
    fy0 = np.zeros(shape=len(x))
    fz0 = np.zeros(shape=len(x))
    N = len(x)
    dt = deltat
    #mass = mass
    
    #update the position at t+dt
    for i in prange(N):
        x[i] = (x[i] + vx[i] * dt + fx[i] * dt * dt * 0.5 / mass ) % (xhi-xlo)
        y[i] = (y[i] + vy[i] * dt + fy[i] * dt * dt * 0.5 / mass ) % (yhi-ylo)
        z[i] = ( z[i] + vz[i] * dt + fz[i] * dt * dt * 0.5 / mass)# % (zhi-zlo)
        
    # save the force at t
    fx0 = fx.copy()         # tested it also only with fx and it still worked
    fy0 = fy.copy()         # tested it also only with fy and it still worked
    fz0 = fz.copy()         # tested it also only with fz and it still worked
    # update acceleration at t+dt
    fx, fy, fz, epot = force.forceLJ(x, y, z, xlo, xhi, ylo, yhi, zlo, zhi, eps, sigma, cutoff, epsWall, cutoffwall, sigWall)

    # update the velocity
    for i in prange(N):        
        vx[i] += 0.5 * dt * (fx[i] + fx0[i]) / mass
        vy[i] += 0.5 * dt * (fy[i] + fy0[i]) / mass
        vz[i] += 0.5 * dt * (fz[i] + fz0[i]) / mass
    
    return x, y, z, vx, vy, vz, fx, fy, fz, epot

 
@njit(parallel=True)
def KineticEnergy(vx, vy, vz, mass):

# calcualte the kinetic energy in joule
    ekin = 0
    N = len(vx)
    i = 0
    
    for i in prange(N):
        ekin += 0.5 * mass * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i]* vz[i])
    return ekin


######### added by Jonas ###########
@njit
def calc_rho(Ngr, hist):
    rho = np.zeros(len(hist))
    for i in range(len(hist)):
        rho[i] = hist[i] / settings.N / Ngr

    return rho


@njit
def update_hists(hist_x, hist_y, hist_z, x,y,z, dr, N):
    for i in prange(N-1):
        # j = i + 1
        for j in range(i+1,N):
            xij = force.pbc(x[i], x[j], settings.xlo, settings.xhi) # calculate pbc distance
            yij = force.pbc(y[i], y[j], settings.ylo, settings.yhi) # can never exceed l/2
            zij = abs(z[i] - z[j])         # can never exceed 2l
 
            if xij != 0:
                bin = int(xij / dr) # find the bin
                hist_x[bin] += 2 # we are counting pairs
            if yij != 0:
                bin = int(yij / dr) # find the bin
                hist_y[bin] += 2 # we are counting pairs
            if zij != 0:
                bin = int(zij / dr) # find the bin
                hist_z[bin] += 2 # we are counting pairs

    return hist_x, hist_y, hist_z