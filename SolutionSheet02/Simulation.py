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
    fx, fy, fz, _,_,epot = force.forceLJ(x, y, z, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi, settings.eps, settings.sig, settings.cutoff,  settings.epsWall, settings.cutoffwall, settings.sigwall)


    # open documents for eq run
    if wT:
        fileoutputeq = open(Tname+ str(everyN)+ 'eq', "w")
        fileoutputEn = open(Tname+ str(everyN)+ 'Energyeq', "w")
        misc.WriteTrajectory3d(fileoutputeq, 0,x,y,z) # noch anpassen an L, Lz
        misc.WriteEnergy(fileoutputEn, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

    
    # do eq run
    for i in tqdm(range(settings.eqsteps)):
        x, y, z, vx, vy, vz, fx, fy, fz, _, _, _ = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi, settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass,  settings.epsWall, settings.cutoffwall, settings.sigwall)
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

    f0sum = 0
    f2lsum = 0
    #open shit for prod
    if wT:
        fileoutputprod = open(Tname+ str(everyN)+ 'prod', "w")
        fileoutputEnprod = open(Tname+ str(everyN)+ 'Energyprod', "w")
        misc.WriteTrajectory3d(fileoutputprod, 0,x,y,z) # noch anpassen an L, Lz
        misc.WriteEnergy(fileoutputEnprod, 0, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))

    # do prod
    Ngr = 0
    nbins_xy = int(settings.l / settings.drxy)
    nbins_z = int(settings.l * 2 / settings.dr)

    hist_x = np.zeros(nbins_xy) 
    hist_y = np.zeros(nbins_xy) 
    hist_z = np.zeros(nbins_z)  # taking as many bins per distance as for x and y
    for i in tqdm(range(1,settings.nsteps)):
        x, y, z, vx, vy, vz, fx, fy, fz, f0, fl, _ = update.VelocityVerlet(x, y, z, vx, vy, vz, fx, fy, fz, settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.zlo, settings.zhi, settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass,  settings.epsWall, settings.cutoffwall, settings.sigwall)
        f0sum += f0
        f2lsum += fl
       # save shit

        if i % everyN == 0:
            if wT:
                misc.WriteTrajectory3d(fileoutputprod, i,x,y,z) # noch anpassen an L, Lz
                misc.WriteEnergy(fileoutputEnprod, i, epot, update.KineticEnergy(vx, vy, vz,settings.mass), np.sum(vx**2), np.sum(vy**2), np.sum(vz**2))
        
        if i % settings.nAnalyze == 0:

            hist_x,hist_y, hist_z = update.update_hists(hist_x,hist_y, hist_z, x, y, z , settings.dr, settings.drxy, settings.N)

            Ngr += 1 # another position
    rho_x = update.calc_rho(Ngr, hist_x, settings.drxy)
    rho_y = update.calc_rho(Ngr, hist_y, settings.drxy)
    rho_z = update.calc_rho(Ngr, hist_z, settings.dr)

    f0 = settings.mass*f0sum /(settings.nsteps-1) /settings.l  /settings.l /settings.N * (-1)# normalize by particles, steps and surface (and surface vector !)
    fl = settings.mass*f2lsum /(settings.nsteps-1) /settings.l  /settings.l  /settings.N

    return rho_x, rho_y, rho_z, f0, fl



