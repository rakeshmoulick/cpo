"""
centre_fluc_E.py
-------------------------------------------------------
OBJ: To observe the fluctuation of electric field at the domain centre  
-------------------------------------------------------
RUN AS: python3 dispersion_npz_EF.py
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

file_name = 'potloc.txt'
path = '../data/files/'
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    time, E = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
    print('No data')
    exit()  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
eps0 = constants('electric constant')
e = constants('elementary charge')
me = constants('electron mass')
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([100,100/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
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
fig,ax = plt.subplots(1,1,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# ----------------------------------------------------------------------------
ax.plot(time, E, 'r',linestyle='-',linewidth=1.0,label='$EF$')
ax.set_xlabel('$\omega_{pe}t$')
ax.set_ylabel('$Electric \; Field$')
ax.set_title('Electric Field at the Domain Centre')
ax.legend(loc='upper right',framealpha=0.5)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.savefig(pjoin(path, 'flucE.jpg'),dpi=1200)

plt.show()
