import numpy as np
import settings
import initialize
import update
import force
import misc
from tqdm import tqdm
from numba import njit, prange

def SimluationSheet4(wT,Tname, everyN):
    # initialize shit
    settings.init()
    x, y, z, vx, vy, vz = initialize.InitializeAtoms()
    fx, fy, fz, epot = force.forceLJ(x, y, z, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi, settings.eps, settings.sig, settings.cutoff,  settings.epsWall, settings.cutoffwall, settings.sigwall)


    # open documents for eq run
    if wT:
        fileoutputeq = open(Tname+ str(everyN)+ 'eq', "w")
        fileoutputEn = open(Tname+ str(everyN)+ 'Energyeq', "w")
        misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
        misc.WriteEnergy(fileoutputEn, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

    
    # do eq run
    for i in tqdm(range(settings.eqsteps)):
        x, y, z, vx, vy, vz, fx, fy, fz, _ = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi, settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass,  settings.epsWall, settings.cutoffwall, settings.sigwall)
        # save shit every n
        if i % everyN == 0:
            if wT:
                misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
                misc.WriteEnergy(fileoutputEn, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

        # if i % 50:
        #     print(f'x, y, z, vx, vy, vz, fx, fy, fz = {x[0], y[0], z[0], vx[0], vy[0], vz[0], fx[0], fy[0], fz[0]}')

        # rescale
        if i % 10 == 0:
            vx, vy, vz = initialize.rescalevelocity(vx, vy, vz,settings.Tdesired, initialize.temperature(vx, vy, vz))

    #open shit for prod
    if wT:
        fileoutputprod = open(Tname+ str(everyN)+ 'prod', "w")
        misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z)


    # do prod
    Ngr = 1
    nbins_xy = int(settings.l/2 / settings.dr)
    nbins_z = int(settings.zhi*2 / settings.dr)

    hist_x = np.zeros(nbins_xy) 
    hist_y = np.zeros(nbins_xy) 
    hist_z = np.zeros(nbins_z)  # taking as many bins per distance as for x and y
    for i in tqdm(range(1,settings.nsteps)):
        x, y, z, vx, vy, vz, fx, fy, fz, _ = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi, settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass,  settings.epsWall, settings.cutoffwall, settings.sigwall)
       # save shit
        if i % everyN == 0:
            if wT:
                misc.WriteTrajectory3d(fileoutputprod, i,x,y,z)
        
        if i % settings.nAnalyze == 0:

            hist_x,hist_y, hist_z = update.update_hists(hist_x,hist_y, hist_z, x, y, z , settings.dr, settings.N)
            rho_x = update.calc_rho(Ngr, hist_x)
            rho_y = update.calc_rho(Ngr, hist_y)
            rho_z = update.calc_rho(Ngr, hist_z)
            Ngr += 1 # another position
    
    return rho_x, rho_y, rho_z



