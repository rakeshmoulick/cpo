"""
Plots the kinetic energy of the total system along with cold, hot and beam electrons.
Appropriate files must be chosen using 'file_name' and 'path'.
Use %reset -f to clear the python workspace.
Data File Invoked: ke_1024.txt
Run as:
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os.path
from os.path import join as pjoin
from scipy.constants import value as constants
import sys

file_name = 'ke.txt'
path = '../data/files/'
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    t, kei, kee = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
    print('No data')
    exit()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
eps0 = constants('electric constant')
e = constants('elementary charge')
me = constants('electron mass')
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ke = kei + kee # Total kinetic energy of the electrons
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
fig,ax = plt.subplots(1,3,figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# vline1 = 38
# vline2 = 300
ax[0].plot(t,ke,'r',linestyle='-',linewidth=1.0,label='KE')
ax[0].set_xlabel('$\omega_{pe}t$')
ax[0].set_ylabel('K_{tot}')
ax[0].legend(loc='upper right',framealpha=0.5)
# ax[0].axvline(vline1,linestyle='--',color='black',linewidth=2.0)
# ax[0].axvline(vline2,linestyle='--',color='black',linewidth=2.0)

ax[1].plot(t,kee,'g',linestyle='-',linewidth=1.0, label='$KE_{C}$')
ax[1].set_xlabel('$\omega_{pe}t$')
ax[1].set_ylabel('$K_{e}$')
ax[1].legend(loc='upper right',framealpha=0.5)
# ax[1].axvline(vline1,linestyle='--',color='black',linewidth=2.0)
# ax[1].axvline(vline2,linestyle='--',color='black',linewidth=2.0)

ax[2].plot(t,kei,'b',linestyle='-',linewidth=1.0, label='$KE_{H}$')
ax[2].set_xlabel('$\omega_{pe}t$')
ax[2].set_ylabel('$K_{i}$')
ax[2].legend(loc='upper right',framealpha=0.5)
# ax[2].axvline(vline1,linestyle='--',color='black',linewidth=2.0)
# ax[2].axvline(vline2,linestyle='--',color='black',linewidth=2.0)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#plt.savefig(pjoin(path,'ke_%d.png'%(vd)),dpi=dpi)

#fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
#ax.loglog(t,ke,'b',linestyle='-',linewidth=1.0)

plt.show()
