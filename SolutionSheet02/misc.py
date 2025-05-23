import settings
import numpy as np
import math
import sys


def WriteEnergy(fileenergy, itime, epot, ekin, vx2, vy2, vz2):
    
    fileenergy.write("%i %e %e %e %e %e\n" % (itime, epot, ekin, vx2, vy2, vz2))

def WriteTrajectory(fileoutput, itime, x, y):

    fileoutput.write("ITEM: TIMESTEP \n")
    fileoutput.write("%i \n" % itime)
    fileoutput.write("ITEM: NUMBER OF ATOMS \n")
    fileoutput.write("%i \n" % (settings.n1*settings.n2))
    fileoutput.write("ITEM: BOX BOUNDS \n")
    fileoutput.write("%e %e xlo xhi \n" % (settings.xlo, settings.xhi))
    fileoutput.write("%e %e xlo xhi \n" % (settings.ylo, settings.yhi))
    fileoutput.write("%e %e xlo xhi \n" % (-1, 1))
    fileoutput.write("ITEM: ATOMS id type x y z \n")
    
    for i in range(0, len(x)):
        z = 0
        fileoutput.write("%i %i %e %e %e \n" % (i, i, (x[i]%settings.xhi), (y[i]%settings.yhi), z))
        

def inputset():
    return settings.xlo, settings.xhi, settings.ylo, settings.yhi, settings.eps, settings.sig, settings.cutoff, settings.deltat, settings.mass
    
def squarevelocity(vx, vy, mass):
    vx2 = 0
    vy2 = 0
    i = 0
    for i in range(0, len(vx)):
        vx2 += vx[i]**2
        vy2 += vy[i]**2
    return 0.5*mass*vx2, 0.5*mass*vy2


def WriteTrajectory3d(fileoutput, itime, x, y, z):

    fileoutput.write("ITEM: TIMESTEP \n")
    fileoutput.write("%i \n" % itime)
    fileoutput.write("ITEM: NUMBER OF ATOMS \n")
    fileoutput.write("%i \n" % (settings.N))
    fileoutput.write("ITEM: BOX BOUNDS \n")
    fileoutput.write("%e %e xlo xhi \n" % (0, settings.l))
    fileoutput.write("%e %e xlo xhi \n" % (0, settings.l))
    fileoutput.write("%e %e xlo xhi \n" % (0, settings.l*2))
    fileoutput.write("ITEM: ATOMS id type x y z \n")
    
    for i in range(0, settings.N):
        fileoutput.write("%i %i %e %e %e \n" % (i, i, x[i] % settings.l, y[i] % settings.l, z[i] % (settings.l*2)))
    
    
