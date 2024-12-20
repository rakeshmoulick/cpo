"""
den_fluc.py
-------------------------------------------------------
OBJ: To plot the density fluctuations
-------------------------------------------------------
RUN AS: python3 den_fluc.py ../input.ini
--------------------------------------------------------
DEVELOPER: DR. RAKESH MOULICK
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import math as mt
import os.path
from os.path import join as pjoin
from scipy.constants import value as constants
import sys

file_name_e = 'denFlucE.txt'
file_name_i = 'denFlucI.txt'
path = '../data/files/'
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name_e)):
    time, dne = np.loadtxt(pjoin(path,file_name_e),unpack=True)
else:
    print('No data')
    exit()
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name_i)):
    time, dni = np.loadtxt(pjoin(path,file_name_i),unpack=True)
else:
    print('No data')
    exit()    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
eps0 = constants('electric constant')
e = constants('elementary charge')
me = constants('electron mass')
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([200,200/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot the figure
fig,ax = plt.subplots(3,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# ----------------------------------------------------------------------------
# vline1 = 38
# vline2 = 300
# ax.axvline(vline1,linestyle='--',color='black',linewidth=2.0)
# ax.axvline(vline2,linestyle='--',color='black',linewidth=2.0)
# ----------------------------------------------------------------------------
ax[0].plot(time, (dne), 'r',linestyle='-',linewidth=1.0,label='$\delta n_{e}$')
ax[0].set_xlabel('$\omega_{pe}t$')
ax[0].set_ylabel('$\delta n_{e}$')
ax[0].legend(loc='upper right',framealpha=0.5)

ax[1].plot(time, (dni), 'g',linestyle='-',linewidth=1.0,label='$\delta n_{i}$')
ax[1].set_xlabel('$\omega_{pe}t$')
ax[1].set_ylabel('$\delta n_{i}$')
ax[1].legend(loc='upper right',framealpha=0.5)

ax[2].plot(time, np.abs(dni-dne), 'b',linestyle='-',linewidth=1.0,label='$\delta  n$')
ax[2].set_xlabel('$\omega_{pe}t$')
ax[2].set_ylabel('$\delta n_{d} = (\delta n_{i}-\delta n_{e})$')
ax[2].legend(loc='upper right',framealpha=0.5)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

plt.savefig(pjoin(path, 'fluctuation.jpg'),dpi=1200)

plt.show()
