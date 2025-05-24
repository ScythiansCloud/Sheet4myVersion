#settings

# velocity: m*s^-1
# position: m
# acceleration: ms^-2
# energy: joule
# temperature: K

import numpy as np

def init():
    global eqsteps
    eqsteps = 50000
    global nsteps            # number of time step to analyze
    nsteps = 100000

    global nAnalyze
    nAnalyze = 3

    global mass              # mass of the LJ particles (gram/mole)
    mass = 39.95
    global kb                # boltzmann's constant (kcal/mole/K) 
    kb =0.0019858775 * 4.184e-6
    global Tdesired          # temperature of the experiment in K
    Tdesired = 300.
    global eps               # eps in LJ (kcal/mole)
    eps = 0.29788162 * 4.184e-6
    global sig              # r0 in LJ (nm)
    sig = 0.188
    global cutoff            # cutoff arbitrary at 2.5 r0
    cutoff = 2.5*sig
    global deltat            # time step (fs)
    deltat = 2
    global epsWall
    epsWall = 1.4887* 4.184e-6
    global cutoffwall
    cutoffwall =0.0376 * 2.5
    global sigwall
    sigwall = 0.0376

    
    # number of particle = n1*n2 distributed on s square lattice
    global n1
    n1 = 6
    global n2
    n2 = 6
    global n3
    n3 = 12
    global N
    N = n1*n2*n3

    # desired density
    global rho
    rho = 0.25 * sig**-3 # / \sigma^2
    global l
    l = (N/(2*rho))**(1/3)

    # for histogramms
    global dr 
    dr = sigwall/30
    global drxy
    drxy = l/100
    

    # box size
    global xlo
    xlo = 0
    global xhi
    xhi = l
    global ylo
    ylo = 0
    global yhi
    yhi = l
    global zlo
    zlo = 0
    global zhi
    zhi = 2*l
    
    global deltaxy # lattice parameter to setup the initial configuration on a lattice
    deltaxy = (xhi - xlo)/n1
    global deltaz
    deltaz = (zhi - zlo)/n3
    
    #rescaling of temperature
    global Trescale
    Trescale = 1 #1 = rescale temperature; 0 = no rescaling
    
    

    
