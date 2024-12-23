"""
plot_results.py
-------------------------------------------------------
OBJ: To plot the electric potential, electric field and 
density fields of electrons and ions 
-------------------------------------------------------
RUN AS: python3 plot_results.py ../input.ini
-------------------------------------------------------
DEVELOPER: DR. RAKESH MOULICK
"""
import numpy as np
from numpy import diff
import matplotlib.pyplot as plt # type: ignore
import matplotlib as mp
from mpl_toolkits import mplot3d
from scipy.constants import value as constants
from scipy.constants import pi
import os.path
from os.path import join as pjoin
from celluloid import Camera # type: ignore
import sys
import configparser 
# ---------------------------------------------
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
eV = e
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
config_path = sys.argv[1]
# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)
# ---------------IMPORT REQUIRED DATA FROM INPUT.INP-----------
NUM_TS = config.getint('time','NUM_TS')
DT_coeff = config.getfloat('diagnostics','DT_coeff')
write_interval = config.getint('diagnostics','write_interval')
n0 = config.getfloat('population','density')
mi = config.getfloat('population','massI')
mi = mi*AMU
NC = config.getint('grid','NC')
nParticlesE = config.getint('population','nParticlesE')
nParticlesI = config.getint('population','nParticlesI')
# ------------------------------------------------------
ni0 = n0 
ne0 = n0

file_name = 'processed_results_all.npz'
path = '../data/files/'

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1
Lx = 2*np.pi
L = Lx/NC
wp = np.sqrt(n0*e**2/(eps0*me)) # Total Plasma Frequency
DT = DT_coeff*(1.0/wp)
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    data = np.load(pjoin(path,file_name))
    x = data['x']
    ndi = data['ndi']
    nde = data['nde']
    rho = data['rho']
    phi = data['phi']   
    EF = data['EF']
else:
    print('No data')
    exit()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = x.reshape(DATA_TS,(NC+1))
xx = x[0,:] # SPATIAL GRID: Select only a column, all the columns have the same value
dx = xx[1]-xx[0] # Normalized StepSize

# create array of the w_pe*t
ts = np.arange(0, NUM_TS+write_interval, write_interval)*DT_coeff

# Reshape Data
# -----------------------------------------------------------------------------
EF = EF.reshape(DATA_TS,(NC+1))
rho = rho.reshape(DATA_TS,(NC+1))
phi = phi.reshape(DATA_TS,(NC+1))
nde = nde.reshape(DATA_TS,(NC+1))
ndi = ndi.reshape(DATA_TS,(NC+1))
# -----------------------------------------------------------------------------
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([300, 60]) #Figure size in mm (FOR SINGLE FIGURE)
figsize1 = np.array([100, 60])
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
k=0
fig, ax = plt.subplots(1,3, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# ----------------------------------------------------------------------------
# camera = Camera(fig)
MaxStp = DATA_TS
while k<MaxStp:
    wpet = k*write_interval*DT_coeff

    ax[0].cla()
    ax[0].plot(x[k,:], phi[k,:], 'k-')
    # ax[0].set_ylim(-2.5, 2.5)
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('$\phi$')    
    ax[0].set_title(k)
    ax[0].grid(True)
    
    ax[1].cla()
    ax[1].plot(x[k,:], EF[k,:],'k-')
    # ax[1].set_ylim(-0.25, 0.25)
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('$E$')
    ax[1].set_title('$\omega_{pe}t = $'+ str(wpet))
    ax[1].grid(True)

    ax[2].cla()
    ax[2].plot(x[k,:], nde[k,:], 'b-', label='$n_{e}$')
    ax[2].plot(x[k,:], ndi[k,:], 'r-', label='$n_{i}$')
    # ax[2].set_ylim(-2, 2)
    ax[2].set_xlabel('x')
    ax[2].set_ylabel('Density')
    #ax[2].set_title(k)
    ax[2].legend(loc='upper right')
    ax[2].grid(True)

    # camera.snap()
    # print(k)
    # print(nde[k,0], nde[k,1], nde[k,-1], nde[k,-2])
    k = k + 20
    plt.pause(0.1)

plt.show()
# ---------------------------------------------------------------------------
# animation = camera.animate()
# animation.save('phi.mp4', fps=5.0)	
# ----------------------------------------------------------------------------
# Path-2: If a continuous running plot is required.
#----------------------------------------------------
# for i in range(DATA_TS):
#     # print(i)
#     #ax[0].cla()
#     ax[0].plot(x[i,:], phi[i,:])
#     ax[0].set_ylim(-1000, 1000)
#     ax[0].set_xlabel('x')
#     ax[0].set_ylabel('$\phi$')
#     #ax[0].set_title(i)
    
#     #ax[1].cla()
#     ax[1].plot(x[i,:], EF[i,:])
#     ax[1].set_ylim(-20, 20)
#     ax[1].set_xlabel('x')
#     ax[1].set_ylabel('$E$')
#     #ax[1].set_title(i)
    
#     #ax[2].cla()
#     ax[2].plot(x[i,:], ndi[i,:], 'r-', label='$n_{i}$')
#     ax[2].plot(x[i,:], ndec[i,:], 'b--', label='$n_{e1}$')
#     ax[2].plot(x[i,:], ndeh[i,:], 'k--', label='$n_{e2}$')
#     ax[2].set_ylim(-1, 3)
#     ax[2].set_xlabel('x')
#     ax[2].set_ylabel('Density')
#     #ax[2].set_title(i)
#     ax[2].legend(loc='upper right')
#     camera.snap()
#     #plt.pause(0.1)

#animation = camera.animate()
#animation.save('phi.mp4', fps=5.0)
# ----------------------------------------------------------------------------
# Path-3: If a static plot is required.
# ---------------------------------------
# indx = 450    
# ax[0].plot(x[indx,:], phi[indx,:])
# ax[0].set_ylim(-1000, 1000)
# ax[0].set_xlabel('x')
# ax[0].set_ylabel('$\phi$')
# #ax[0].set_title(indx)

# ax[1].plot(x[indx,:], EF[indx,:])
# ax[1].set_ylim(-20, 20)
# ax[1].set_xlabel('x')
# ax[1].set_ylabel('$E$')
# #ax[1].set_title(indx)

# ax[2].plot(x[indx,:], ndi[indx,:], 'r-', label='$n_{i}$')
# ax[2].plot(x[indx,:], ndec[indx,:], 'b--', label='$n_{e1}$')
# ax[2].plot(x[indx,:], ndeh[indx,:], 'k--', label='$n_{e2}$')
# ax[2].set_xlabel('x')
# ax[2].set_ylabel('$n_{i}$, $n_{e}$')
# ax[2].legend(loc=0)
# plt.savefig(pjoin(path,'sheath.png'),dpi=dpi)

# ----------------------------------------------------------------------------
# plot: plot the density contours
# ----------------------------------------------------------------------------
# print(ts.shape)
# print(xx.shape)
# print(EF.shape)
# fig, ax = plt.subplots(figsize=figsize1/25.4,constrained_layout=True,dpi=ppi)
# Z = phi
# plt.pcolor(xx, ts, Z, cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))
# cbar = plt.colorbar()
# cbar.set_label('$\phi$')
# ax.set_xlabel('x')
# ax.set_ylabel('$\omega_{pe}t$')
# plt.savefig(pjoin(path,'phiContour.png'),dpi=dpi)
# ------------------------------------------------------
# fig, ax = plt.subplots(figsize=figsize1/25.4,constrained_layout=True,dpi=ppi)
# Z = ndec
# # plt.pcolor(xx, ts, Z, cmap='viridis',shading='auto',vmin=np.min(Z),vmax=np.max(Z))
# plt.pcolor(xx, ts, Z, cmap='viridis',shading='auto',vmin=0,vmax=1)
# cbar = plt.colorbar()
# cbar.set_label('$n_{e}$')
# ax.set_xlabel('x')
# ax.set_ylabel('$\omega_{pe}t$')
# plt.savefig(pjoin(path,'denContour.png'),dpi=dpi)

# ------------------------------------------------------
# Z = nde
# pos, time = np.meshgrid(xx, ts)
# fig = plt.figure(figsize=figsize1/10,constrained_layout=True,dpi=ppi)
# ax = plt.axes(projection='3d')
# ax.plot_surface(pos, time, Z, cmap=plt.cm.cividis, edgecolor='none')
# ax.view_init(elev=80, azim=40)
# ax.set_xlabel('x')
# ax.set_ylabel('$\omega_{pe}t$')
# #ax.set_zlabel('$n_{e}$')
# ax.zaxis.set_tick_params(labelbottom=False)
# ax.set_zlim(0,1.5)
# plt.savefig(pjoin(path,'denSurf.png'))
