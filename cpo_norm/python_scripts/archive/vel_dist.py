"""
# This file plots the particle data at fixed times. This is the file to generate plots for the paper.
# Only citing the path to the data file is enough to plot,
# No specific input corresponding to parameteric variation is required
# Run as: 
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import pandas as pd
import seaborn as sns
from scipy.constants import value as constants 
from scipy.interpolate import make_interp_spline

sns.set(style='whitegrid')
import os.path
from os.path import join as pjoin
import sys

#path = sys.argv[1]
path = '../data/'
path_fig = './'
# ------------- Comments -------------------------------------------------------
# Path to data file for vd=20 is
# ../data/particle_data/part_vd_20/data
# Path to data file for vd=80 is
# ../data/particle_data/part_vd_80/data
# Figures to be saved in ../data/particle_data/part_vd_20/figs for vd=20 &
# ../data/particle_data/part_vd_80/figs for vd=80
#-------------------------------------------------------------------------------
# Calculate the equilibrium densities of the beam, cold and hot electrons
alp = 1.0
n0 = 1E10
ni0 = n0
nec0 = n0/(1+alp)
neh0 = alp*nec0

# print(neb0)
# print(nec0)
# print(neh0)
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS = 50000
write_interval = 100
DT_coeff = 0.01
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
n0 = 1E10
Te = 1*e   # cold electron temperature in joule
Ti = 0.026*e
mi = 40*AMU
# ----------------------------------------------------------------------------
NC = 512
dx = 1
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
LDH = np.sqrt(eps0*Te/(neh0*e**2)) # Hot electron Debye Length
LDC = np.sqrt(eps0*Te/(nec0*e**2)) # Cold electron Debye Length
LD = np.sqrt(eps0*Te/(n0*e**2)) # Characteristic Debye length

wpi = np.sqrt(ni0*e**2/(eps0*mi))
wp  = np.sqrt(n0*e**2/(eps0*me)) # Total Plasma Frequency

wpec = np.sqrt((nec0*e**2)/(eps0*me))
wpeh = np.sqrt(neh0*e**2/(eps0*me))

nParticlesE = 20000
nParticlesI = nParticlesE 
xl = NC*dx

electron_cold_spwt = (nec0*xl*LD)/(nParticlesE)
electron_hot_spwt = (neh0*xl*LD)/(nParticlesE)
# print(electron_cold_spwt, electron_hot_spwt)
# Define the time at which data is required
DT_coeff = 0.01
#write_interval = 200
# File index value is k(say), change this index to get the plot at required time.
# Calculate the wpet for a k-value by using DT_coeff beforehand.
k = [0, 45000]
nbins = 100
count = np.zeros([len(k), nbins])
vele = np.zeros([len(k), nbins])
wpet = np.zeros(len(k))

for i in range(len(k)):   
    wpet[i] = k[i]*DT_coeff # Value of wpet
    
    file_name_ec = 'ec%d.txt'%(int(k[i]))
    file_name_eh = 'eh%d.txt'%(int(k[i]))
    
    # Load Cold electron data
    if os.path.exists(pjoin(path,file_name_ec)):
        xec,vec = np.loadtxt(pjoin(path,file_name_ec),unpack=True)
    else:
        print('No data')
        exit()
    
    # Load Hot electron data
    if os.path.exists(pjoin(path,file_name_eh)):
        xeh,veh = np.loadtxt(pjoin(path,file_name_eh),unpack=True)
    else:
        print('No data')
        exit()
    
    v = np.concatenate((vec, veh)) 
    
    mini = np.min(v)-0.5
    maxi = np.max(v)+0.5
    #nbins = 51
    delta = (maxi - mini)/(nbins)
    
    bin = np.zeros(nbins)
    for j in range(len(v)):
        c = int(np.floor((v[j] - mini)/delta))
        if j<20000:    
            bin[c] += 1*electron_cold_spwt
        else:
            bin[c] += 1*electron_hot_spwt

    y = np.linspace(np.min(v), np.max(v), len(bin))
    
    
    count[i,:] = bin
    vele[i,:] = y
    
    
# spline0 = make_interp_spline(vele[0,:], count[0,:])
# spline1 = make_interp_spline(vele[1,:], count[1,:])

# X_0 = np.linspace(vele[0,:].min(), vele[0,:].max(), 500)
# X_1 = np.linspace(vele[1,:].min(), vele[1,:].max(), 500)

# Y_0 = spline0(X_0)
# Y_1 = spline1(X_1)

# plt.plot(X_0, Y_0/np.max(Y_0), 'r-', linewidth=2, label='$\u03C9_{pe}t=%d$'%(wpet[0]))
# plt.plot(X_1, Y_1/np.max(Y_1), 'k--', linewidth=2, label='$\u03C9_{pe}t=%d$'%(wpet[1]))

maxc = np.max(count[0,:])
plt.plot(vele[0,:], count[0,:]/maxc, 'r--',linewidth=2, label='$\u03C9_{pe}t=%d$'%(wpet[0]))
plt.plot(vele[1,:], count[1,:]/maxc, 'k-', linewidth=2, label='$\u03C9_{pe}t=%d$'%(wpet[1]))

plt.xlabel('v')
plt.ylabel('f(v)')
plt.legend(loc=0)
#plt.xlim(-30, 30)
#plt.ylim(0, 1)

# Save Figure
plt.savefig(pjoin(path_fig,'veldist.png'),format='png',dpi=600)

plt.show()
